# Mitochondrial variant calling pipeline

Reference mitochondrial genome used:

https://www.ncbi.nlm.nih.gov/nuccore/CM008198.1

The mitochondrial variant calling pipeline is used to the determine which SNPs are present in each of the input samples mitochondrion. The pipeline consists of the following steps:

1. Preliminary quality control check of the data.
2. Adapter trimming of the reads.
3. Quality trimming of the reads.
4. Second quality control check, to see effects of trimming.
5. Aligning the reads to the reference genome.
6. Mark duplicates.
7. Variant calling.

the end result from the pipeline is a collection of vcf files containing the SNP's found on the mitochondrial genome.

## 1. Preliminary quality control check of the data.

The first step in the pipeline is to do some quality control(QC) checks of the submitted fastq files. For this purpose FastQC[put reference here] is used. 

## 2. Adapter trimming of the reads.

The second step the in the pipeline is to trim the adapters of the reads, for this purpose fastp[put reference here] is used. Fastp automatically detects the adapters and trims them therefore only needing the fastq files as input.

## 3. Quality trimming of the reads.

After adapter trimming another round of trimming is done, but based on quality of the reads. This is done using trimmomatic[reference here], a java based tool used for illumina data. The following settings are used:

- MINLEN:35
- LEADING:20
- TRAILING:20 
- SLIDINGWINDOW:3:15 
- AVGQUAL:20

## 4. Second quality control check, to see effects of trimming.

After quality trimming another round of QC checks is done to see the effects of the pre-proccessing on the raw data. Once again fastqc is used.

## 5. Aligning the reads to the reference genome.

With the quality trimming done the preprocessing steps have been finished and the reads can now be mapped against a reference genome(see top of the readme). The tool used for this step is BWA[reference here] and in paticular BWA-MEM.  

## 6. Mark duplicates.

In this step duplicate reads are marked as they do not hold any relevant additional information and can introduce a bias in the variant calling step. To mark the duplicates picard[reference here] is used. 

## 7. Variant calling.

With duplicates marked the next step is to call the SNP's in the sample. This is done with GATK Haplotypecaller[reference here].

To run the pipeline the following command can be used, after changing the yaml file to match your correct input, intermediar, and output folders, to check whether everything is setup correctly.

	snakemake -n -s {name of the snakefile}


If that doens't throw an error, you can run the pipeline as shown below.

	snakemake -c {amount of cores} -s {name of the snakefile}

The yaml file is specified at the top of the snakefile and contains the following keys that have to be set correctly, with paths going from your work directory:

- fastq_folder: Folder containing the raw fast files.
- fastqc_folder: Folder where the first round of QC reports will be stored.
- fastqc_folder2: Folder where the second round of QC reports will be stored.
- ref: The refrence genome to be used (full path).
- aling: The folder in which the alignment will be stored.
- vcf: The folder in which the vcf files will be stored.