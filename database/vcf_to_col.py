#Script to turn the vcf downloaded from EVA into the right format for the tables in the database

import pandas as pd
import numpy as np

df = pd.read_csv("ChrMT-Run8.txt", sep='\t', header=0)
df_breed = pd.read_csv("../retrieve_data_pipeline/result/test_info_EVA.csv", sep='\t', header=0)

ref = list(df['REF'])
alt = list(df['ALT'])
pos = list(df['POS'])
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
with open('tables/accession_table.txt', 'w') as f:
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
with open('tables/snp_table.txt', 'w') as f:
	for num, i in enumerate(snp_col):
		f.write('%s|%s|%s|%s|%s\n' % (num + 1, bs_id[num], ref_id[num], i, hetro_homo_alt[num]))

#refsite table
with open('tables/refsite_table.txt', 'w') as f:
	for num, i in enumerate(ref):
		f.write('%s|%s|%s\n' % (int(num + 1), int(pos[num]), str(i)))


		