import sqlite3 as sq
from sqlite3 import Error
import argparse
import pandas as pd
import math

parser = argparse.ArgumentParser(prog = 'db_stuff.py', description = 'does db stuff.')
parser.add_argument('-i', '--input', help='path to database.', required=True )
parser.add_argument('-v', '--vcf', help='path to vcf file.', required=True )  
parser.add_argument('-msa', '--msa', help='path to MSA output file.', required=True )
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

query = "SELECT s.accession_id, a.biosample, a.breed, a.country, count(s.accession_id) FROM snps AS s INNER JOIN refsites as r ON s.refsite_id = r.refsite_id INNER JOIN accessions AS a ON s.accession_id = a.accession_id WHERE r.position = %s AND  s.alt_allele = '%s' " % (pos_snp[0], snp_col[0])
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

for col in df_all.columns:
	one_snp = ''
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
				one_snp = one_snp + ref + '/'
				if ',' in alt:
					split_alt = alt.split(',')
					one_snp = one_snp + split_alt[(int(split_gt[1]) -1)]
				else:
					one_snp = one_snp + alt

			else:
				print(alt, gt)
				split_alt = alt.split(',')
				one_snp = split_alt[(int(split_gt[0]) -1)] + '/' + split_alt[(int(split_gt[1]) -1)]	

		else:
			one_snp = one_snp + alt + '/' + alt
		
		parchement_snp.append(one_snp)
		col_list.append(df_all[col])
		parchement_ref.append('%s/%s' % (ref,ref))


df = pd.concat(col_list, axis=1)
df.columns = df.iloc[0]
df = df.drop(df.index[0])
df.iloc[0] = parchement_ref
print(len(parchement_ref), len(parchement_snp))

# for i in accession ids
bs = list(zip(*acc_info))
bs_list = []
br_list = []
ct_list = []
print('extract correct snps')
for num, i in enumerate(bs[0]):
	#	get the snps
	#	fill a row in the df with the snps and where there is no snp put ref
	row = []


	try:
		cursor.execute('SELECT refsite_id, alt_allele, hom_het FROM snps WHERE accession_id = %i;' % (i))
		snp_info = cursor.fetchall()
		snp = list(zip(*snp_info))
		for column in df.columns:
			one_snp = ''
			if column in snp[0]:
				gt = str(snp[2][snp[0].index(column)])
				split_gt = gt.split('/')

				ref = df[column].iloc[0][0]
				alt = snp[1][snp[0].index(column)]

				if split_gt[0] != split_gt[1]:
					if split_gt[0] == '0':
						one_snp =  ref + '/' + alt
					else:
						print(alt, gt, column)
						split_alt = alt.split(',')
						one_snp = split_alt[0] + '/' + split_alt[1]
				else:
					one_snp =  alt + '/' + alt
			
			else:
				#this was '.', changed it to show the ref nt df[column][0]
				one_snp =  df[column].iloc[0][0] + '/' + df[column].iloc[0][0]

			row.append(one_snp)
		df.loc[len(df)] = row
		bs_list.append(bs[1][num])
		br_list.append(bs[2][num])
		ct_list.append(bs[3][num])


	except Error as e:
		print(f"The error '{e}' occurred")


drop = []
print(df)

df_drop = df.iloc[0:len(df)]
for num, column in enumerate(df_drop.columns):
	#print(df_drop[column].unique())
	#when doing dots and '.' in df_drop[column].unique()
	if len(df_drop[column].unique()) < 2:	
		drop.append([num,column])
	elif int(parchement_dp[num]) < args.dp:
		drop.append([num,column])
	elif float(parchement_qual[num]) < args.qual:
		drop.append([num,column])
	elif len(parchement_snp[num]) != 3:
		drop.append([num,column])
	else:
		one = False
		for i in df_drop[column].unique():
			if len(i) != 3:
				one=True
		if one == True:
			drop.append([num,column])

retract = 0

for i in drop:	
	df = df.drop(i[1], axis=1)
	del parchement_snp[i[0]-retract]
	del parchement_dp[i[0]-retract]
	del parchement_qd[i[0]-retract]
	del parchement_qual[i[0]-retract]
	retract +=1

print(df)

for i in range(len(parchement_snp)):
	print('%s\t%s\t%s\t%s' % (parchement_snp[i], parchement_qual[i],parchement_dp[i], parchement_qd[i]))


#create biosample column

bs_list.insert(0, 'Reference cow')
df['biosample'] = bs_list
#create breed column

br_list.insert(0, 'Reference cow')
df['breed'] = br_list
#create country column

ct_list.insert(0, 'UK')
df['country'] = ct_list

#need to do something to get the name there
parchement_snp.append('Parchment  ')
parchement_snp.append('Parchment  ')
parchement_snp.append('Parchment  ')
df.loc[len(df)] = parchement_snp


#set biosample as index
#df = df.set_index('biosample')
print('going to csv')
#df = df.iloc[:, : 10]
df.to_csv(args.msa, sep='\t')

