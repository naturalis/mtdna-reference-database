import pandas as pd
import numpy as np
import argparse

parser = argparse.ArgumentParser(prog = 'create_query.py', description = 'This script creates a query that counts the amount of matching snps the sample has in common with all samples in the database.')
parser.add_argument('-i', '--input', help='The input vcf file for the script.', required=True )           

args = parser.parse_args()
df = pd.read_csv(args.input, sep='\t', names =['#CHROM', 'POS', 'REF', 'ALT', 'QUAL', 'SAMPLE'])

ref = list(df['REF'])
alt = list(df['ALT'])
pos = list(df['POS'])
biosamples = list(df.columns)[5:]

snp_col = []
pos_snp = []


for num, i in enumerate(biosamples):
	one_sample = list(df[i])
	for j in range(len(one_sample)):
		split_sp = one_sample[j].split(':')

		if split_sp[0] in ['0/0', './.']:
			continue
		else:
			gt = split_sp[0]
			split_gt = gt.split('/') 

			if int(split_gt[0]) > 0 and split_gt[0] != split_gt[1]:
				alt_nsps = alt[j].split(',')
				snp_col.append('%s,%s' % (alt_nsps[int(split_gt[0])-1], alt_nsps[int(split_gt[1])-1]))
			else:
				snp_col.append(alt[j].split(',')[int(split_gt[1])-1])

			pos_snp.append(pos[j])



with open('query.txt', 'w') as f:
	string = "SELECT accession_id, count(*) FROM snps AS s INNER JOIN refsites as r ON s.refsite_id = r.refsite_id WHERE r.position = %s AND  s.alt_allele = '%s' " % (pos_snp[0], snp_col[0])
	#create WHERE ... AND statement that checks for every position matching with the right
	#snp and counts the amount matching for every accesion
	for i, snp in enumerate(snp_col[1:]):
		string = string + "OR r.position = %s AND s.alt_allele = '%s' " % (pos_snp[i+1], snp)
		if i + 1 == len(snp_col[1:]):
			string = string + "OR r.position = %s AND s.alt_allele = '%s' GROUP BY accession_id ORDER BY count();" % (pos_snp[i+1], snp)
		
	f.write(string)
		