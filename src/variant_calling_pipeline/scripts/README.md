# Scripts

### get_range.py

The get_range.py script takes a total of two commandline arguments:

    -i  --input 
                      Bed file containing all the regions in which the al.
    -o  --output 
                      The output file to write the ranges to.

The script takes the bed file and loops over every line keeping track of the positions, if there is an interuption in the continuety of the positions a new region is defined. This keeps going until the end of the file is reached.

### vcf_filter.py

The vcf_filter.py script takes a total of three commandline arguments:

    -i  --input 
                      The path to the input VCF file.
    -o  --output 
                      The path and common prefix for to the output files, in which the results will be stored
    -r  --regions
    				  The path to the regions bed file.

The script loops through the SNPs in the VCF file and checks for multiple criteria on which it filters the SNPs into different output files. The criteria are on whether or not it's whithin the region with a coverage of 10, QUAL score, DP value, and GQ score.