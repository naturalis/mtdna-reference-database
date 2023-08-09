import pysam
import logging
import argparse
import sys
import re
import os
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment


def parse_bed(bed_file_location):
    """
    Parse a BED file and return a list of intervals as dictionaries.
    :param bed_file_location: Path to the BED file
    :return: List[dict]: List of intervals in the format [{"chrom": contig, "start": pos, "end": pos}, ...]
    """
    logging.info(f"Going to read BED file from {bed_file_location}")
    intervals = []
    with open(bed_file_location, 'r') as bedfile:
        for line in bedfile:
            if not line.startswith("#"):  # skip comments if any
                parts = line.strip().split('\t')
                chrom = parts[0]
                start = int(parts[1])
                end = int(parts[2])
                intervals.append({"chrom": chrom, "start": start, "end": end})
    logging.debug(intervals)
    return intervals


def parse_breeds_table(file_location):
    """
    Parses the TSV file that contains the curated list of breeds. Returns the IDs of the breeds, i.e. the first
    column sans header.
    :param file_location: Path to the TSV file
    :return: Dictionary: List of sample IDs
    """
    logging.info(f"Going to read breeds table from {file_location}")
    breeds = {}
    with open(file_location, 'r') as file:
        next(file)  # Skip header
        for line in file:
            stripped_line = line.strip()  # Remove leading/trailing whitespaces
            if stripped_line:  # Check if line is not blank
                fields = stripped_line.split('\t')  # Get fields
                breeds[fields[0]] = fields[2]
    logging.debug(breeds)
    return breeds


def parse_parchment_vcf(vcf_file, sequences):
    logging.info(f"Going to read parchment VCF from {vcf_file}")

    # "constant" indices in the VCF
    pos_idx = 1
    ref_idx = 3
    alt_idx = 4

    # generate sample name and initialize list
    basename = os.path.basename(vcf_file)  # Gets "P1_vcf.vcf.gz"
    stem = os.path.splitext(basename)[0]  # Gets "P1_vcf.vcf"
    sample = os.path.splitext(stem)[0]  # Gets "P1_vcf"
    sequences[sample] = []

    # start reading file
    with open(vcf_file, 'r') as file:
        for line in file:

            # skip header lines
            if line.startswith('#'):
                continue

            # skip blank lines
            stripped_line = line.strip()
            if stripped_line:
                fields = stripped_line.split('\t')
                allele = fields[alt_idx]

                # grow the list to wherever we're going to insert
                pos = int(fields[pos_idx])
                while len(sequences[sample]) <= pos:
                    sequences[sample].append('N')

                match = re.match(r'^(.+?),', allele)
                if match:
                    result = match.group(1)
                    sequences[sample][pos] = result
                else:
                    sequences[sample][pos] = fields[ref_idx]

    return sequences

def interpolate_sequences(vcf_file, reference_file, samples=None, contig='M', sequences=None):
    """
    Interpolates the single reference sequence with variants from the VCF file. Optionally does this only on a subset of
    samples from the VCF file. Optionally can add the result to an existing sequences dictionary, in order to merge
    files.
    :param vcf_file: Path to the bgzipped, tabixed VCF file
    :param reference_file: Path to the faidxed FASTA file
    :param samples: List of samples to only look at (default uses all)
    :param contig: Name of contig to look at (default is 'M', i.e. how the mitome is called in this data set)
    :param sequences: Dictionary where key is sample, value is list of nucleotides
    :return: Dictionary: Keys are sample IDs, values are lists of sites.
    """
    # Load reference sequence and VCF
    with pysam.FastaFile(reference_file) as ref:
        logging.info(f"Going to read reference sequence from {reference_file}")
        reference_sequence = ref.fetch(contig)
    logging.info(f"Going to read SNPs from VCF file {vcf_file}")
    vcf = pysam.VariantFile(vcf_file)

    # Create samples list if not defined. This is variable so that in one pass we can use the curated list
    # of breeds, and in another pass we let the script extract the ID of the parchment sample.
    if samples is None:
        samples = list(vcf.header.samples)
        logging.debug(f"No list of samples provided, read from VCF header: {samples}")
    else:
        logging.debug(f"Provided list of samples: {samples}")

    # Dictionary to store sequences for each sample. This is variable so that in one pass we instantiate the
    # dictionary, and in the next we grow it by adding the parchment sample.
    if sequences is None:
        logging.info("Initializing new set of sequences")
        sequences = {sample: list(reference_sequence) for sample in samples}
    else:
        logging.info("Expanding existing set of sequences")

    # Iterate over variants (SNPs etc.)
    for record in vcf:

        # Iterate over samples for that SNP, possibly restricting ourselves to the curated set.
        for sample, sample_data in record.samples.items():
            if sample in samples:
                indices = sample_data.allele_indices
                if len(indices) > 0:
                    index = indices[0]  # It's mtDNA so we treat it as homozygous.
                    if index is not None:
                        sequences[sample][record.pos] = record.alleles[index]

    return sequences


def pad_indels(sequences):
    """
    Attempts to deal with indels by padding the shorter variants with '-'. This will likely give shitty
    alignments, but maybe we are lucky and there are no multi-allelic indels (famous last words).
    :param sequences: Dictionary where key is sample, value is list of nucleotides
    :return: Padded sequence dictionary
    """

    # calculate size of matrix
    nchar = None
    for sample in sequences.keys():
        if nchar is None:
            nchar = len(sequences[sample])
            continue
        if nchar != len(sequences[sample]):
            nn = len(sequences[sample])
            logging.error(f"{nchar} != {nn}")

    for i in range(nchar):

        # calculate the length of the longest variant at this site
        maxlen = 1
        for sample in sequences.keys():
            if len(sequences[sample][i]) > maxlen:
                maxlen = len(sequences[sample][i])

        # pad others if they are shorter than maxlen
        if maxlen > 1:
            for sample in sequences.keys():
                if len(sequences[sample][i]) < maxlen:
                    padding = maxlen - len(sequences[sample][i])
                    gap = '-' * padding
                    sequences[sample][i] += gap
    return sequences


def filter_sites(sequences, intervals):
    """
    Reduces the interpolated sequences to only those sites that overlap with the provided intervals.
    :param sequences: Dictionary where key is sample, value is list of nucleotides
    :param intervals: Intervals from the BED file
    :return: Filtered sequence dictionary
    """
    # make an empty copy, to which we will append the sites inside the intervals
    filtered = {sample: [] for sample in sequences.keys()}
    for interval in intervals:
        start = interval["start"]
        end = interval["end"]  # BED end coordinate is exclusive, but so is python list slicing!!!
        for sample in sequences.keys():
            subseq = sequences[sample][start:end]
            filtered[sample].extend(subseq)
    return filtered


def write_output(sequences, outformat):
    """
    Writes the result to stdout in one of the formats that AlignIO understands, e.g. phylip or nexus
    :param sequences: Dictionary where key is sample, value is list of nucleotides
    :param outformat: String that specifies the AlignIO output, e.g. phylip or nexus
    """
    alignment = MultipleSeqAlignment([])
    for sample in sequences.keys():
        seq = ''.join(sequences[sample])
        sr = SeqRecord(seq, id=sample)
        sr.annotations["molecule_type"] = "DNA"
        alignment.append(sr)
    AlignIO.write(alignment, sys.stdout, outformat)


def main():
    parser = argparse.ArgumentParser(description="Process BED and VCF files.")

    # add arguments
    parser.add_argument("--bed", required=True, help="Location of the BED intervals file.")
    parser.add_argument("--snp-vcf", required=True, help="Location of the VCF SNP panel file.")
    parser.add_argument("--parchment-vcf", required=True, help="Location of the VCF parchment sample file.")
    parser.add_argument("--breeds", required=True, help="Location of breeds table.")
    parser.add_argument("--ref", required=True, help="Location of reference FASTA.")
    parser.add_argument("--outformat", required=True, default="nexus", help="Output format.")
    parser.add_argument("--verbose", action="store_true", help="Increase logging verbosity.")
    args = parser.parse_args()

    # configure logging
    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    # prepare return values
    bed_file = args.bed
    vcf_snp_file = args.snp_vcf
    vcf_parchment_file = args.parchment_vcf
    breeds_file = args.breeds
    ref_file = args.ref
    outformat = args.outformat

    # log results
    logging.debug(f"BED file location: {bed_file}")
    logging.debug(f"VCF SNP panel file location: {vcf_snp_file}")
    logging.debug(f"VCF parchment sample file location: {vcf_parchment_file}")
    logging.debug(f"FASTA reference sequence location: {ref_file}")
    logging.debug(f"TSV breeds file location: {breeds_file}")
    logging.debug(f"Output format: {outformat}")

    return bed_file, vcf_snp_file, vcf_parchment_file, breeds_file, ref_file, outformat


if __name__ == "__main__":
    # preprocess input files from command line arguments
    bed_f, snp_f, parchment_f, breeds_f, ref_f, outformat = main()
    breeds = parse_breeds_table(breeds_f)
    intervals = parse_bed(bed_f)

    # call interpolation, then merge the parchment into it using homegrown VCF parser :-(
    sequences = interpolate_sequences(snp_f, ref_f, list(breeds.keys()))
    sequences = parse_parchment_vcf(parchment_f, sequences)

    # filter the interpolated sequences to retain intervals
    sequences = filter_sites(sequences, intervals)

    # pad the indels
    #sequences = pad_indels(sequences)

    # write output
    #write_output(sequences, outformat)
    for sample in sequences.keys():
        print(f">{sample}")
        seq = ''.join(sequences[sample])
        print(f"{seq}")
