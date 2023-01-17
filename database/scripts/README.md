# Scripts

### vcf_to_col.py

The vcf_to_col.py script takes five command line inputs:


    -i   --input INPUT
                          The input vcf file for the script.
    -t1  --table1 TABLE1
                          The output file to write the accessions table to.
    -t2  --table2 TABLE2
                          The output file to write the snps table to.
    -t3  --table3 TABLE3
                          The output file to write the refrence site table to.
    -m   --metadata METADATA
                          The info file containing metadata on the cows.


The script takes a vcf file and the metadata file created by the cow_metadata_pipeline as its -i and -m inputs respectively and uses the information within these files to create the tables for the database.

### create_nexus_splitstree.py

The create_nexus_splitstree.py script takes the following commandline arguments:

    -i   --input 
                path to database. 
    -v   --vcf 
                path to vcf file.  
    -msa --msa 
                path to the MSA output file.  
    -nex --nexus 
                path to the NEXUS output file.  
    -d   --dp 
                DP threshold to filter on.    
    -q   --qual 
                QUAL threshold to filter on.   


The script takes the database and the VCF file of the parchment filtered on the regions with a depth of 10 as it's input. The script then querries the database for matches; filters the matches on the QUAL and DP values; and then formats the results in the NEXUS fromat. The resulting NEXUS file can then be used as input for splitstree

### two_seq_matrix.py

The two_seq_matrix.py script takes the following commandline arguments:

    -i   --input 
                path to database. 
    -v   --vcf 
                path to vcf file.  
    -msa --msa 
                path to the MSA output file.  
    -f   --fasta 
                path to the fasta output file.  
    -d   --dp 
                DP threshold to filter on.    
    -q   --qual 
                QUAL threshold to filter on.   

The script is almost a copy of the create_nexus_splitstree.py script, the only difference is that the an MSA with two sequences per cow is created and from these two sequences the one containing all the SNPs from the cow is saved in a fasta file. This fasta file is used as input for the haplotype_network.R script.

### for_r_msa.py

The for_r_msa.py script takes the following commandline arguments:

    -i   --input 
                path to database. 
    -v   --vcf 
                path to vcf file.  
    -msa --msa 
                path to the MSA output file.    
    -d   --dp 
                DP threshold to filter on.    
    -q   --qual 
                QUAL threshold to filter on.   

The script is also almost a copy of the create_nexus_splitstree.py script. The only difference once again lies in the way the MSA looks. The resulting CSV file is also used as input for the haplotype_network.R script.

### haplotype_network.R

As the name already implies this R script is used to create the haplotype networks. It takes the results from the for_r_msa and two_seq_matrix python scripts as input.