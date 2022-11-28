import pandas as pd
import numpy as np
import argparse

parser = argparse.ArgumentParser(prog = 'create_query.py', description = 'This script creates a query that counts the amount of matching snps the sample has in common with all samples in the database.')
parser.add_argument('-i', '--input', help='The input vcf file for the script.', required=True )           

args = parser.parse_args()
df = pd.read_csv(args.input, sep='\t', header=0)

ref = list(df['REF'])
alt = list(df['ALT'])
pos = list(df['POS'])
biosamples = list(df.columns)[9:]

snp_col = []
ref_snp = []
pos_snp = []


for num, i in enumerate(biosamples):
	one_sample = list(df[i])
	for j in range(len(one_sample)):
		split_sp = one_sample[j].split(':')

		if split_sp[0] in ['0']:
			continue
		else:
			split_gt = split_sp[0]

			if int(split_gt) > 1:
				split_alt = alt[j].split(',')
				split_rat = np.array(split_sp[1].split(',')[1:]) 
				snp_col.append(split_alt[np.argmax(split_rat)])

			else:
				snp_col.append(alt[j])

			ref_snp.append(ref[j])
			pos_snp.append(pos[j])


with open('query.txt', 'w') as f:
	string = "SELECT accession_id, count(*) FROM snps AS s INNER JOIN refsites as r ON s.refsite_id = r.refsite_id WHERE r.position = %s AND  s.alt_allele = '%s' " % (pos_snp[0], snp_col[0])
	#create WHERE ... AND statement that checks for every position matching with the right
	#snp and counts the amount matching for every accesion
	for i, snp in enumerate(snp_col[1:]):
		string = string + "OR r.position = %s AND s.alt_allele = '%s' " % (pos_snp[i+1], snp)
		if i + 1 == len(snp_col[1:]):
			print('here')
			string = string + "OR r.position = %s AND s.alt_allele = '%s' GROUP BY accession_id ORDER BY count();" % (pos_snp[i+1], snp)
		
	f.write(string)
print(string)
		