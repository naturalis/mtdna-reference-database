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

The script takes a vcf file and the metadata file created by the cow_metadata_pipeline as its -i and -m inputs respectively and uses the information within these files to create the tables for the database.

### create_query.py

The create_query.py script takes one command line input:

'''
  -i INPUT, --input INPUT
                        The input vcf file for the script.
'''

The script takes a vcf file of a parchment sample as input and then creates a query that will count the amount of matching snps each cow in the database has with the sample.