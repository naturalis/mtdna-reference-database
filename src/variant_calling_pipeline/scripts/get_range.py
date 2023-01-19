#get the positions in which the depth of a bam file is over 10 into a txt file 		samtools depth in.bam | awk '$3 > 5'
#loop over the lines
#if the  position[i] is equal to position[i - 1]
# then keep going
#otherwise position[i] is the new start position and position [i-1] is the end position
#that can then be added to a list

import pandas as pd 
import argparse

parser = argparse.ArgumentParser(prog = 'Get_ranges.py', description = 'This script retrieves the ranges to be extracted from the bam file, that have a depth of 10 or higher.')
parser.add_argument('-i', '--input', help='The input file for the script. Containing the positions and depth of 10 or higher in the bam file (bed file).', required=True )           
parser.add_argument('-o', '--output', help='The output file for the script.', required=True )           
args = parser.parse_args()


df = pd.read_csv(args.input, sep='\t', names=['chr', 'pos', 'dp'])

start = df['pos'][0]
end = 0
ranges = ''
#codes = {'NC_037328.1':'chr1', 'NC_037329.1':'chr2', 'NC_037330.1':'chr3', 'NC_037331.1':'chr4', 'NC_037332.1':'chr5', 'NC_037333.1':'chr6', 'NC_037334.1':'chr7', 'NC_037335.1':'chr8', 'NC_037336.1':'chr9', 'NC_037337.1':'chr10', 'NC_037338.1':'chr11', 'NC_037339.1':'chr12', 'NC_037340.1':'chr13', 'NC_037341.1':'chr14', 'NC_037342.1':'chr15', 'NC_037343.1':'chr16', 'NC_037344.1':'chr17', 'NC_037345.1':'chr18', 'NC_037346.1':'chr19', 'NC_037347.1':'chr20', 'NC_037348.1':'chr21', 'NC_037349.1':'chr22', 'NC_037350.1':'chr23', 'NC_037351.1':'chr24', 'NC_037352.1':'chr25', 'NC_037353.1':'chr26', 'NC_037354.1':'chr27', 'NC_037355.1':'chr28', 'NC_037356.1':'chr29', 'NC_037357.1':'chrX'}

for num, i in enumerate(df['pos']):
	if i - 1 == df['pos'].iloc[num - 1]:
		continue

	elif i - 1 != df['pos'].iloc[num - 1]:
		if df['chr'].iloc[num - 1] == 'NC_006853.1':
			continue
		else:
			end = df['pos'].iloc[num - 1]
			chrom = df['chr'].iloc[num - 1]
			ranges = ranges + "%s\t%s\t%s\n" % (chrom, start, end)
			start = i

with open(args.output, 'w') as f:
	f.write(ranges)

