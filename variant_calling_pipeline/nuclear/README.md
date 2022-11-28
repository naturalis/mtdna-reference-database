# Nuclear variant calling pipeline

Refrence nucleair genome used:

https://www.ncbi.nlm.nih.gov/assembly/GCF_002263795.1/

The nuclear variant calling pipeline is used to the determine which SNPs are present is each of the input samples chromosome regions that have a coverage of higher then ten. The pipeline consists of the following steps:

1. Preliminary quality control check of the data.
2. Adapter trimming of the reads.
3. Quality trimming of the reads.
4. Aligning the reads to the reference genome.
5. Second quality control check, to see effects of trimming.
6. Retrieving reagions with a depth of at least ten.
7. Mark duplicates.
8. Variant calling.
9. Vcf filtering.

the end result from the pipeline is a collection of text files containing filtered vcf records based on different criteria. Such as, SNPs with a DP (depth) of at least 10, a quality score of at least 20, a GQ score of at least 20, all of the above and SNPs that are within the regions the alignment says the depth is at least 10.

## 1. Preliminary quality control check of the data.

The first step in the pipeline is to do some quality control checks of the submitted fastq files. For this purpose FastQC[put reference here] is used. 

## 2. Adapter trimming of the reads.

The second step the in the pipeline is to trim the adapters of the reads, for this purpose fastp[put reference here] is used. Fastp automatically detects the adapters and trims them therefore only the 

## 3. Quality trimming of the reads.

## 4. Aligning the reads to the reference genome.

## 5. Second quality control check, to see effects of trimming.

## 6. Retrieving reagions with a depth of at least ten.

## 7. Mark duplicates.

## 8. Variant calling.

## 9. Vcf filtering.


To run the pipeline the following command can be used, after changing the yaml file to match your correct input, intermediar, and output folders, to check whether everything is setup correctly.

'''
snakemake -n -s {name of the snakefile}
''' 

If that doens't throw an error, you can run the pipeline as shown below.

'''
snakemake -c {amount of cores} -s {name of the snakefile}
''' 

The yaml file is specified at the top of the snakefile and contains the following keys:

	- fastq_folder:
	- fastqc_folder:
	- fastqc_folder2:
	- ref:
	- aling:
	- vcf:
	- scripts:
	- dp: