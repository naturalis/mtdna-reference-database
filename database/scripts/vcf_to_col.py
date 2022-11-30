#Script to turn the vcf downloaded from EVA into the right format for the tables in the database

import pandas as pd
import numpy as np
import argparse

parser = argparse.ArgumentParser(prog = 'vcf_to_col.py', description = 'This script extracts information from the vcf file and metadata file to then put it in the correct format for the database tables.')
parser.add_argument('-i', '--input', help='The input vcf file for the script.', required=True )           
parser.add_argument('-t1', '--table1', help='The output file to write the accessions table to.', required=True )
parser.add_argument('-t2', '--table2', help='The output file to write the snps table to.', required=True )
parser.add_argument('-t3', '--table3', help='The output file to write the refrence site table to.', required=True )
parser.add_argument('-m', '--metadata', help='The info file containing metadata on the cows.', required=True )           

args = parser.parse_args()

df = pd.read_csv(args.input , sep='\t', header=0)
df_breed = pd.read_csv(args.metadata, sep='\t', header=0)

ref = list(df['REF'])
alt = list(df['ALT'])
pos = list(df['POS'])
chrom = list(df['#CHROM'])
biosamples = list(df.columns)[9:]

snp_col = []
ref_snp = []
pos_snp = []
hetro_homo_alt = []
bs_id = []
ref_id = []


#Make a list with geographical differnt cows to test 
	#[SAMN05803884, SAMN02142985, SAMN01915355, SAMN13389915, SAMN05788490, SAMEA4644757, SAMN02842704, SAMEA33668668, SAMN02842707, SAMEA5159763, SAMEA33001168, SAMN02843065, SAMN02225744, SAMN06282409, SAMN05199762, SAMN05217649, SAMN09379816]
	#[IR, SC, FR (SW), FRxBE, UK, SW, FR (C), IT, DU, FR (CS), FR (NW), FR (NW), KR, PL, FR (C), FR (E), BE]

#accessions table
with open(args.table1, 'w') as f:
	for num, i in enumerate(biosamples):
		one_sample = list(df[i])
		for j in range(len(one_sample)):
			split_sp = one_sample[j].split(':')

			if split_sp[0] in ['0/0','./.']:
				continue
			else:
				split_gt = split_sp[0].split('/')

				if split_gt[0] != split_gt[1]:
					try:
						split_alt = alt[j].split(',')
						split_rat = np.array(split_sp[1].split(','))
						idx = np.argmax(split_rat)
						if idx == 0:
							continue
						else: 
							snp_col.append(split_alt[idx])
					except:
						snp_col.append(alt[j])

				elif split_gt[0] == split_gt[1]:
					try:
						split_alt = alt[j].split(',')
						split_rat = np.array(split_sp[1].split(',')) 
						idx = np.argmax(split_rat)
						if idx == 0:
							continue
						else: 
							snp_col.append(split_alt[idx])
					except:
						snp_col.append(alt[j])

				ref_snp.append(ref[j])
				pos_snp.append(pos[j])
				hetro_homo_alt.append(split_sp[0])
				bs_id.append(num + 1)
				ref_id.append(j + 1)

		try:
			f.write('%s|%s|%s|%s|%s\n' % (num + 1, i, df_breed['breed'].loc[df_breed['biosample'] == i].iloc[0],
				df_breed['country'].loc[df_breed['biosample'] == i].iloc[0],
				df_breed['direction'].loc[df_breed['biosample'] == i].iloc[0]))
		except:
			continue

#snp table
with open(args.table2, 'w') as f:
	for num, i in enumerate(snp_col):
		f.write('%s|%s|%s|%s|%s\n' % (num + 1, bs_id[num], ref_id[num], i, hetro_homo_alt[num]))

#refsite table
with open(args.table3, 'w') as f:
	for num, i in enumerate(ref):
		f.write('%s|%s|%s|%s\n' % (int(num + 1), chrom[num], int(pos[num]), str(i)))


		