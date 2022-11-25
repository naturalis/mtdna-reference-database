import pandas as pd 
import numpy as np
import itertools
import argparse

parser = argparse.ArgumentParser(prog = 'Get_ranges.py', description = 'This script retrieves the ranges to be extracted from the bam file, that have a depth of 10 or higher.')
parser.add_argument('-i', '--input', help='The input vcf file for the script.', required=True )           
parser.add_argument('-o', '--output', help='The output file to write the results to.', required=True )
parser.add_argument('-r', '--regions', help='The regions.bed file containing the regions where depth was 10 or higher.', required=True )           

args = parser.parse_args()

#to do:
#make sure snps are in defined regions

#codes = {'NC_037328.1':'chr1', 'NC_037329.1':'chr2', 'NC_037330.1':'chr3', 'NC_037331.1':'chr4', 'NC_037332.1':'chr5', 'NC_037333.1':'chr6', 'NC_037334.1':'chr7', 'NC_037335.1':'chr8', 'NC_037336.1':'chr9', 'NC_037337.1':'chr10', 'NC_037338.1':'chr11', 'NC_037339.1':'chr12', 'NC_037340.1':'chr13', 'NC_037341.1':'chr14', 'NC_037342.1':'chr15', 'NC_037343.1':'chr16', 'NC_037344.1':'chr17', 'NC_037345.1':'chr18', 'NC_037346.1':'chr19', 'NC_037347.1':'chr20', 'NC_037348.1':'chr21', 'NC_037349.1':'chr22', 'NC_037350.1':'chr23', 'NC_037351.1':'chr24', 'NC_037352.1':'chr25', 'NC_037353.1':'chr26', 'NC_037354.1':'chr27', 'NC_037355.1':'chr28', 'NC_037356.1':'chr29', 'NC_037357.1':'chrX'}

df_regions = pd.read_csv(args.regions, sep='\t', names=['chr', 'start', 'end'])

chrom_bool = False
length = []
 
with open(args.input, 'r') as f_in:
	with open(args.output + "_dp10.txt", 'w') as f_out:
		with open(args.output + "_qual20.txt", 'w') as f_out2:
			with open(args.output + "_gq10.txt", 'w') as f_out3:
				with open(args.output + "_all.txt", 'w') as f_out4:
					with open(args.output + "_in_regions.txt", 'w') as f_out5:
						for line in f_in:

							if line[:6] == "#CHROM":
								chrom_bool=True
								continue

							if chrom_bool:
								split_line = line.split('\t')
								df2 = df_regions.loc[(df_regions['chr'] == split_line[0]) & (df_regions['start'] < int(split_line[1])) & (df_regions['end'] > int(split_line[1]))]
								if len(df2) != 0:
									f_out5.write("%s\t%s\t%s\t%s\t%s\t%s" % (split_line[0], split_line[1], split_line[3], split_line[4], split_line[5], split_line[-1]))
								try: 
									if float(split_line[5]) >= 19.99:
										f_out2.write("%s\t%s\t%s\t%s\t%s\t%s" % (split_line[0], split_line[1], split_line[3], split_line[4], split_line[5], split_line[-1]))
									if int(split_line[-1].split(':')[2]) > 9:
										f_out.write("%s\t%s\t%s\t%s\t%s\t%s" % (split_line[0], split_line[1], split_line[3], split_line[4], split_line[5], split_line[-1]))	
									if int(split_line[-1].split(':')[3]) > 9:
										f_out3.write("%s\t%s\t%s\t%s\t%s\t%s" % (split_line[0], split_line[1], split_line[3], split_line[4], split_line[5], split_line[-1]))
									if float(split_line[5]) >= 19.99 and int(split_line[-1].split(':')[2]) > 9 and int(split_line[-1].split(':')[3]) > 9:
										length.append(split_line[0])
										f_out4.write("%s\t%s\t%s\t%s\t%s\t%s" % (split_line[0], split_line[1], split_line[3], split_line[4], split_line[5], split_line[-1]))
								except IndexError:
									print("IndexError: SNP cannot be called GT field contains \"%s\"" % split_line[-1]) 
							else:
								continue

print(len(length))



			

