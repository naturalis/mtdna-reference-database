import sqlite3 as sq
from sqlite3 import Error
import argparse
import pandas as pd
import math

parser = argparse.ArgumentParser(prog = 'db_stuff.py', description = 'does db stuff.')
parser.add_argument('-i', '--input', help='path to database.', required=True )
parser.add_argument('-v', '--vcf', help='path to vcf file.', required=True )  
parser.add_argument('-msa', '--msa', help='path to full MSA output file.', required=True )
parser.add_argument('-f', '--fasta', help='path to output fasta file.', required=True )
parser.add_argument('-d', '--dp', help='DP threshold to filter on.', required=True )  
parser.add_argument('-q', '--qual', help='QUAL threshold to filter on.', required=True )  


args = parser.parse_args()

df_vcf = pd.read_csv(args.vcf, sep='\t',  names =['#CHROM', 'POS', 'REF', 'ALT', 'QUAL', 'SAMPLE'])


ref = list(df_vcf['REF'])
alt = list(df_vcf['ALT'])
pos = list(df_vcf['POS'])
biosamples = list(df_vcf.columns)[5:]

snp_col = []
pos_snp = []

print('create query')
#create query
for num, i in enumerate(biosamples):
	one_sample = list(df_vcf[i])
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

query = "SELECT s.accession_id, a.biosample, a.breed , count(s.accession_id) FROM snps AS s INNER JOIN refsites as r ON s.refsite_id = r.refsite_id INNER JOIN accessions AS a ON s.accession_id = a.accession_id WHERE r.position = %s AND  s.alt_allele = '%s' " % (pos_snp[0], snp_col[0])
#create WHERE ... AND statement that checks for every position matching with the right
#snp and counts the amount matching for every accesion
for i, snp in enumerate(snp_col[1:]):
	query = query + "OR r.position = %s AND s.alt_allele = '%s' " % (pos_snp[i+1], snp)
	if i + 1 == len(snp_col[1:]):
		query = query + "OR r.position = %s AND s.alt_allele = '%s' GROUP BY s.accession_id ORDER BY count(s.accession_id) DESC LIMIT 30;" % (pos_snp[i+1], snp)


#connect to database
try:
    connection = sq.connect(args.input)
    print("Connection to SQLite DB successful")
except Error as e:
    print(f"The error '{e}' occurred")


cursor = connection.cursor()

# get the acccessioon_id, biosamples and save in a list, can be index
try:
    cursor.execute('%s' % query)
    acc_info = cursor.fetchall()

except Error as e:
    print(f"The error '{e}' occurred")

#get the reference table
try:
    cursor.execute('SELECT position, refsite_id, refallele FROM refsites;')
    ref_nts = cursor.fetchall()

except Error as e:
    print(f"The error '{e}' occurred")

#make a df with every positon in the reference having its own column
#make the first row the reference 
ref = list(zip(*ref_nts))
df_all = pd.DataFrame(columns=ref[0])
df_all.loc[-1] = ref[2]
df_all.loc[-2] = ref[1]  
df_all.index = df_all.index + 1  
df_all = df_all.sort_index()

#get the columns in which the parchement has a snp
#and save the snps to add to the dataframe later
col_list = []
parchement_snp = []
parchement_ref = []
parchement_dp = []
parchement_qd = []
parchement_qual = []
second_parchment = []

#split the SNPs into two seperate ones, for hetrozygous SNPs all the mutations are
#put on the same sequence. 
for col in df_all.columns:
	if col in list(df_vcf['POS']):
		p_sp = df_vcf.loc[df_vcf['POS'] == col]['SAMPLE'].iloc[0]
		split_sp = p_sp.split(':')
		gt = split_sp[0]
		parchement_dp.append(split_sp[2])
		parchement_qd.append(split_sp[3])
		split_gt = gt.split('/')

		ref = df_vcf.loc[df_vcf['POS'] == col]['REF'].iloc[0]
		alt = df_vcf.loc[df_vcf['POS'] == col]['ALT'].iloc[0]
		qual = df_vcf.loc[df_vcf['POS'] == col]['QUAL'].iloc[0]

		parchement_qual.append(qual)

		if split_gt[0] != split_gt[1]:
			if split_gt[0] == '0':
				second_parchment.append(ref)
				if int(split_gt[1]) > 1:
					split_alt = alt.split(',')
					parchement_snp.append(split_alt[(int(split_gt[0]) -1)])
				else:
					parchement_snp.append(alt)

			else:
				print(alt, gt)
				split_alt = alt.split(',')
				second_parchment.append(split_alt[(int(split_gt[0]) -1)])
				parchement_snp.append(split_alt[(int(split_gt[1]) -1)])	

		else:
			second_parchment.append(alt)
			parchement_snp.append(alt)
		
		col_list.append(df_all[col])
		parchement_ref.append(ref)

print(len(parchement_snp), len(second_parchment))

df = pd.concat(col_list, axis=1)
df.columns = df.iloc[0]
df = df.drop(df.index[0])
print(len(parchement_ref), len(parchement_snp))
for i in range(len(df.columns)):
	#print(i, len(parchement_snp), len(df.columns))
	if len(df[df.columns[i]].iloc[0]) > len(parchement_snp[i]):
		p_ref = parchement_snp[i]
		p_snp = parchement_snp[i]
		new_p_snp = p_snp + df[df.columns[i]].iloc[0][len(p_ref):]
		parchement_snp[i] = new_p_snp

#split the SNPs into two seperate ones, for hetrozygous SNPs all the mutations are
#put on the same sequence. 
bs = list(zip(*acc_info))
bs_list = []
br_list = []
print('extract correct snps')
for num, i in enumerate(bs[0]):
	#	get the snps
	#	fill a row in the df with the snps and where there is no snp put ref
	row = []
	row2 = []
	try:
		cursor.execute('SELECT refsite_id, alt_allele, hom_het FROM snps WHERE accession_id = %i;' % (i))
		snp_info = cursor.fetchall()
		snp = list(zip(*snp_info))
		for column in df.columns:
			if column in snp[0]:
				gt = str(snp[2][snp[0].index(column)])
				split_gt = gt.split('/')

				ref = df[column].iloc[0]
				alt = snp[1][snp[0].index(column)]

				if split_gt[0] != split_gt[1]:
					if split_gt[0] == '0':
						row2.append(ref)
						row.append(alt)
					else:
						print(alt, gt, column)
						split_alt = alt.split(',')
						row2.append(split_alt[0])
						row.append(split_alt[1])
				else:
					row.append(alt)
					row2.append(alt)
			
			else:
				#this was '.', changed it to show the ref nt df[column][0]
				row.append(df[column][0])
				row2.append(df[column][0])
		df.loc[len(df)] = row2
		df.loc[len(df)] = row
		bs_list.append(bs[1][num])
		bs_list.append(bs[1][num])
		br_list.append(bs[2][num])
		br_list.append(bs[2][num])

	except Error as e:
		print(f"The error '{e}' occurred")


drop = []
print(df)
#filter columns(SNPs) based on whether it's all the same, DP value, and QUAL.
#Very important that the
df_drop = df.iloc[0:len(df)]
for num, column in enumerate(df_drop.columns):
	if len(df_drop[column].unique()) < 2:		
		drop.append([num,column])
	elif int(parchement_dp[num]) < int(args.dp):
		drop.append([num,column])
	elif float(parchement_qual[num]) < int(args.qual):
		drop.append([num,column])
	else:
		one = False
		for i in df_drop[column].unique():
			if len(i) != 1:
				one=True
		if one == True:
			drop.append([num,column])

retract = 0
for i in drop:	
	df = df.drop(i[1], axis=1)
	del parchement_snp[i[0]-retract]
	del second_parchment[i[0]-retract]
	del parchement_dp[i[0]-retract]
	del parchement_qd[i[0]-retract]
	del parchement_qual[i[0]-retract]
	retract +=1

print(df)


#create biosample column
bs_list.insert(0, 'Reference cow')
df['biosample'] = bs_list

#create breed column
br_list.insert(0, 'Reference cow')
df['breed'] = br_list

#add parchement rows
parchement_snp.append('Parchment  ')
parchement_snp.append('Parchment  ')
second_parchment.append('Parchment  ')
second_parchment.append('Parchment  ')
df.loc[len(df)] = second_parchment
df.loc[len(df)] = parchement_snp


print('create MSA')
#make the MSA
for i in df.columns[:-2]:
	ref_len = len(df[i][0])
	idx = 0
	while idx < len(df[i]):

		if ',' in df[i][idx]:
			split_snp = df[i][idx].split(',')

			for num, snp in enumerate(split_snp):

				row = list(df.loc[idx])
				ind = list(df.columns).index(i)
				row[ind] = snp

				if num == 0:
					df.loc[idx] = row
				else:
					df.loc[len(df)] = row
		idx += 1
	for j in range(len(df[i])):
		max_len = df[i].str.len().max()
		if len(df[i][j]) == max_len and df[i][j] != '*':
			idx += 1
			continue
		else:
			diffrence = max_len - len(df[i][j])

			if df[i][j] == '*':
				df[i][j] = '-' + (diffrence * '-')

			else:
				df[i][j] = df[i][j] + (diffrence * '-')
	

#set biosample as index
#df = df.set_index('biosample')
print('going to csv')

df.to_csv(args.msa, sep='\t')
print(df.shape)
with open(args.fasta, 'w') as f:
	for i in range(0,len(df),2):
		row = ''.join(list(df.iloc[i][:-2]))
		f.write('>[%s]\n%s\n' % (i, row))

