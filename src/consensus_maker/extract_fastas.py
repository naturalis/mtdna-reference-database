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


def parse_samples_table(file_location):
    """
    Parses the TSV file that contains the curated list of samples in ../../doc/Breed_Country_Region.tsv. Returns the
    IDs of the samples, i.e. the first column sans header as well as the third column, which contains the broad
    geographic origin ('Superregion').
    :param file_location: Path to the TSV file.
    :return: Dictionary: List of sample IDs.
    """
    logging.info(f"Going to read breeds table from {file_location}")
    samples = {}
    with open(file_location, 'r') as file:
        next(file)  # Skip header
        for line in file:
            stripped_line = line.strip()  # Remove leading/trailing whitespaces
            if stripped_line:  # Check if line is not blank
                fields = stripped_line.split('\t')  # Get fields
                samples[fields[0]] = fields[2]
    logging.debug(samples)
    return samples


def parse_parchment_vcf(vcf_file, sequences):
    """
    Hand-rolled parser of the (single-sample) VCF files coming out of Galaxy. Reads the positions, ref alleles and alt
    alleles to construct an array of nucleotides that is added to the sequences dictionary keyed on the file sten of
    the sample's VCF file.
    :param vcf_file: Location of a single-sample VCF file from Galaxy.
    :param sequences: A dictionary keyed on sample names, whose values are lists of nucleotides.
    :return:
    """
    logging.info(f"Going to read parchment VCF from {vcf_file}")

    # "constant" indices in the VCF
    pos_idx = 1 # genomic startcoordinate
    ref_idx = 3 # reference allele
    alt_idx = 4 # alternative alleles, will assume homozygous because mtDNA

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
                pos = int(fields[pos_idx]) - 1
                alt_allele = fields[alt_idx]
                ref_allele = fields[ref_idx]

                # grow the list to wherever we're going to insert
                while len(sequences[sample]) <= pos:
                    sequences[sample].append('N')

                match = re.match(r'^(.+?),', alt_allele) # e.g. 'G,<*>'
                if match:
                    result = match.group(1)
                    sequences[sample][pos] = result
                else:
                    sequences[sample][pos] = ref_allele
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


def filter_indels(sequences):
    """
    Removes all sites where at least one variant's length != 1 (i.e. indels). This then ought to mean that the result
    does not need to be realigned, and we get around the weird semantics of indels in haplotype networks.
    :param sequences: Dictionary where key is sample, value is list of nucleotides
    :return: Filtered sequence dictionary
    """
    # make an empty copy, to which we will append the sites inside the intervals
    filtered = {sample: [] for sample in sequences.keys()}

    # calculate nchar and do sanity check
    nchar = None
    for sample in sequences.keys():
        if nchar is None:
            nchar = len(sequences[sample])
            continue
        else:
            length = len(sequences[sample])
            if nchar != length:
                logging.error(f"{nchar} != {length} in sample {sample}")

    # iterate over sites while toggling indel boolean to continue to next site
    indel = False
    for i in range(nchar):
        for sample in sequences.keys():
            if len(sequences[sample][i]) != 1:
                indel = True
                break
        if indel:
            indel = False
            continue
        else:
            for sample in sequences.keys():
                base = sequences[sample][i]
                if base in [ 'A', 'C', 'G', 'T']:
                    filtered[sample].append(base)
                else:
                    filtered[sample].append('N')

    return filtered


def main():
    parser = argparse.ArgumentParser(description="Process BED and VCF files.")

    # add arguments
    parser.add_argument("--bed", required=True, help="Location of the BED intervals file.")
    parser.add_argument("--snp-vcf", required=True, help="Location of the VCF SNP panel file.")
    parser.add_argument("--parchment-vcf", required=True, help="Location of the VCF parchment sample file.")
    parser.add_argument("--breeds", required=True, help="Location of breeds table.")
    parser.add_argument("--ref", required=True, help="Location of reference FASTA.")
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

    # log results
    logging.debug(f"BED file location: {bed_file}")
    logging.debug(f"VCF SNP panel file location: {vcf_snp_file}")
    logging.debug(f"VCF parchment sample file location: {vcf_parchment_file}")
    logging.debug(f"FASTA reference sequence location: {ref_file}")
    logging.debug(f"TSV breeds file location: {breeds_file}")

    return bed_file, vcf_snp_file, vcf_parchment_file, breeds_file, ref_file


if __name__ == "__main__":
    # preprocess input files from command line arguments
    bed_f, snp_f, parchment_f, breeds_f, ref_f = main()
    breeds = parse_samples_table(breeds_f)
    intervals = parse_bed(bed_f)

    # call interpolation, then merge the parchment into it using homegrown VCF parser :-(
    sequences = interpolate_sequences(snp_f, ref_f, list(breeds.keys()))
    sequences = parse_parchment_vcf(parchment_f, sequences)

    # filter the interpolated sequences to retain intervals, remove indels
    sequences = filter_sites(sequences, intervals)
    sequences = filter_indels(sequences)

    # write output
    for sample in sequences.keys():
        print(f">{sample}")
        seq = ''.join(sequences[sample])
        print(f"{seq}")
