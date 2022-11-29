# Database

The database folder contains the databases which were created during the project, scripts and a test input. The folder contains the following five folders: 

- data:
		Containing a test vcf file to create a database with. 
- mitochondrial:
		Containing the mitochondrial snp database and its schema.
- nuclear:
		Containing the nuclear snp database (or a link for it) and its schema.
- results:
		Containing the tables created from the information in the vcf file and metadata file, that was used to fill the databases.
- scripts:
		Containing the scripts for creating the data tables that will fill the database and creating a query to retrieve the amout of matched snps. 

# Scripts


### vcf_to_col.py

The vcf_to_col.py script takes five command line inputs:

'''
  -i INPUT, --input INPUT
                        The input vcf file for the script.
  -t1 TABLE1, --table1 TABLE1
                        The output file to write the accessions table to.
  -t2 TABLE2, --table2 TABLE2
                        The output file to write the snps table to.
  -t3 TABLE3, --table3 TABLE3
                        The output file to write the refrence site table to.
  -m METADATA, --metadata METADATA
                        The info file containing metadata on the cows.
'''

The script takes a vcf file and the metadata file created by the cow_metadata_pipeline as its -i and -m inputs and uses the information within these files to create the tables for the database.

### create_query.py

The create_query.py script takes one command line input:

'''
  -i INPUT, --input INPUT
                        The input vcf file for the script.
'''

The script takes a vcf file of a parchment sample as input and then creates a query that will count the amount of matching snps each cow in the database has with the sample.