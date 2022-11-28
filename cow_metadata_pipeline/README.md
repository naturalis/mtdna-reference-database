# Cow_metadata_pipeline

The cow_metadata_pipeline is used to retrieve the breed of the cow samples in either a SRA run file or ENA filereport. The pipeline consists of the following steps:

1. Extracting the biosample IDs from either SRA run file or ENA filereport.
2. Retrieve the discriptions of the biosamples from NCBI database based on the common start of biosample IDs.
3. Extract the actual breed from the retrieved descriptions.
4. Format all the information in a easily accessable manner.

The end result from the pipeline will be a file, for which the name can be specified in the accompanying .yaml file, that contains the breed, biosample ID, experiment ID, run ID, country of origin of the breed, and cardinal direction describing where within the country the breed is from. The information is structured as followed:

|breed|biosample|experiment_id|run|country|direction|
| --- | --- | --- | --- | --- | --- |
|Holstein|SAMN10598568|SRS4148739|SRR8324586|US|WHOLE|
Hereford|SAMN10598569|SRS4148752|SRR8324572|GB|CW|
Vechur|SAMN10131258|SRS3931706|?|?|?|
Belgian Blue|SAMN09379796|SRS3434671|ftp.sra.ebi.ac.uk/vol1/ERZ173/ERZ1738264/ChrMT-Run8-TAUIND-public.vcf.gz,ftp.sra.ebi.ac.uk/vol1/ERZ173/ERZ1738264/Chr2-Run8-TAUIND-public.vcf.gz|BE|WHOLE|
Limousin|SAMN09379797|SRS3434715|ftp.sra.ebi.ac.uk/vol1/ERZ173/ERZ1738264/ChrMT-Run8-TAUIND-public.vcf.gz,ftp.sra.ebi.ac.uk/vol1/ERZ173/ERZ1738264/Chr2-Run8-TAUIND-public.vcf.gz|FR|CS|

To run the pipeline the following command can be used, after changing the yaml file to match your correct input, intermediar, and output folders, to check whether everything is setup correctly.

'''
snakemake -n -s {name of the snakefile}
''' 

If that doens't throw an error, you can run the pipeline as shown below.

'''
snakemake -c 1 -s {name of the snakefile}
''' 

The yaml file is specified at the top of the snakefile and contains the following keys:

	- sra_folder: 
		The folder containing the input SRA run files and/or ENA filereport. ENA files do have to be renamed SraRunTable{name}, as can be seen in the test/input folder.
	- result_folder:
		The folder in which results from the pipeline will be stored. 
	- intermediar_folder:
		The folder in which intermediary files will be stored.
	- breed_folder:
		The folder in which the files containing the breeds from your samples will be stored.
	- breed_file:
		The file containing the different breeds of cow in the used dataset and where they are from. **IMPORTANT: don't delete the file. The file, all_breed_counts.txt in the repository, is manually curated as the NCBI descriptions do not contain country information.** There is also a copy in the backup folder.
	- description_folder:
		The folder in which the raw descriptions from NCBI will be stored.
	- script_folder:
		The folder containing the scripts used in the pipeline.
	- output_file: 
		The name of the output file containing the end result from the pipeline.