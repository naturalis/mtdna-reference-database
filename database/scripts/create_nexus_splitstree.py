#probably for splits tree
import sqlite3 as sq
from sqlite3 import Error
import argparse
import pandas as pd
import math

parser = argparse.ArgumentParser(prog = 'db_stuff.py', description = 'does db stuff.')
parser.add_argument('-i', '--input', help='path to database.', required=True )
parser.add_argument('-v', '--vcf', help='path to vcf file.', required=True )  
parser.add_argument('-msa', '--msa', help='path to output file.', required=True )
parser.add_argument('-nex', '--nexus', help='path to output file.', required=True )
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

#Turn the hetrozygous snps from the parchement VCF files into the correct nucleotide
#to retain the information.
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
				if ref == 'A' and alt == 'G' or ref == 'G' and alt == 'A':
					parchement_snp.append('R')
				elif ref == 'A' and alt == 'C' or ref == 'C' and alt == 'A':
					parchement_snp.append('M') 
				elif ref == 'A' and alt == 'T' or ref == 'T' and alt == 'A':
					parchement_snp.append('W')
				elif ref == 'G' and alt == 'T' or ref == 'T' and alt == 'G':
					parchement_snp.append('K')
				elif ref == 'G' and alt == 'C' or ref == 'C' and alt == 'G':
					parchement_snp.append('S')
				elif ref == 'C' and alt == 'T' or ref == 'T' and alt == 'C':
					parchement_snp.append('Y')
				else:
					parchement_snp.append(alt)
			else:
				split_alt = alt.split(',')
				if split_alt[(int(split_gt[0]) -1)] == 'A' and split_alt[(int(split_gt[1]) -1)] == 'G' or split_alt[(int(split_gt[0]) -1)] == 'G' and split_alt[(int(split_gt[1]) -1)] == 'A':
					parchement_snp.append('R')
				elif split_alt[int(split_gt[0]) -1] == 'A' and split_alt[int(split_gt[1]) -1] == 'C' or split_alt[int(split_gt[0]) -1] == 'C' and split_alt[int(split_gt[1]) -1] == 'A':
					parchement_snp.append('M') 
				elif split_alt[int(split_gt[0]) -1] == 'A' and split_alt[int(split_gt[1]) -1] == 'T' or split_alt[int(split_gt[0]) -1] == 'T' and split_alt[int(split_gt[1]) -1] == 'A':
					parchement_snp.append('W')
				elif split_alt[int(split_gt[0]) -1] == 'G' and split_alt[int(split_gt[1]) -1] == 'T' or split_alt[int(split_gt[0]) -1] == 'T' and split_alt[int(split_gt[1]) -1] == 'G':
					parchement_snp.append('K')
				elif split_alt[int(split_gt[0]) -1] == 'G' and split_alt[int(split_gt[1]) -1] == 'C' or split_alt[int(split_gt[0]) -1] == 'C' and split_alt[int(split_gt[1]) -1] == 'G':
					parchement_snp.append('S')
				elif split_alt[int(split_gt[0]) -1] == 'C' and split_alt[int(split_gt[1]) -1] == 'T' or split_alt[int(split_gt[0]) -1] == 'T' and split_alt[int(split_gt[1]) -1] == 'C':
					parchement_snp.append('Y')
				else:
					parchement_snp.append(alt)
		else:
			parchement_snp.append(alt)
		
		col_list.append(df_all[col])
		parchement_ref.append(ref)


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

# for i in accession ids
bs = list(zip(*acc_info))

print('extract correct snps')
#Turn the hetrozygous snps from the reference database into the correct nucleotide
#to retain the information.
for i in bs[0]:

	#	fill a row in the df with the snps and where there is no snp put ref
	row = []
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
						if ref == 'A' and alt == 'G' or ref == 'G' and alt == 'A':
							row.append('R')
						elif ref == 'A' and alt == 'C' or ref == 'C' and alt == 'A':
							row.append('M') 
						elif ref == 'A' and alt == 'T' or ref == 'T' and alt == 'A':
							row.append('W')
						elif ref == 'G' and alt == 'T' or ref == 'T' and alt == 'G':
							row.append('K')
						elif ref == 'G' and alt == 'C' or ref == 'C' and alt == 'G':
							row.append('S')
						elif ref == 'C' and alt == 'T' or ref == 'T' and alt == 'C':
							row.append('Y')
						else:
							row.append(alt)
					else:
						print(alt, gt, column)
						split_alt = alt.split(',')
						if split_alt[0] == 'A' and split_alt[1] == 'G' or split_alt[0] == 'G' and split_alt[1] == 'A':
							row.append('R')
						elif split_alt[0] == 'A' and split_alt[1] == 'C' or split_alt[0] == 'C' and split_alt[1] == 'A':
							row.append('M') 
						elif split_alt[0] == 'A' and split_alt[1] == 'T' or split_alt[0] == 'T' and split_alt[1] == 'A':
							row.append('W')
						elif split_alt[0] == 'G' and split_alt[1] == 'T' or split_alt[0] == 'T' and split_alt[1] == 'G':
							row.append('K')
						elif split_alt[0] == 'G' and split_alt[1] == 'C' or split_alt[0] == 'C' and split_alt[1] == 'G':
							row.append('S')
						elif split_alt[0] == 'C' and split_alt[1] == 'T' or split_alt[0] == 'T' and split_alt[1] == 'C':
							row.append('Y')
						else:
							row.append(alt)
				else:
					row.append(alt)
			
			else:
				row.append(df[column][0])
		df.loc[len(df)] = row

	except Error as e:
		print(f"The error '{e}' occurred")

drop = []
print(df)

counter_len = 0
counter_uniq = 0
counter_dp = 0
counter_q = 0
df_drop = df.iloc[0:len(df)]
#Filter the snps on different critera, such as the DP and QUAL
for num, column in enumerate(df_drop.columns):
	if len(df_drop[column].unique()) < 2:		
		drop.append([num,column])
		counter_uniq +=1
	elif int(parchement_dp[num]) < args.dp:
		drop.append([num,column])
		counter_dp+=1
	elif float(parchement_qual[num]) < args.qual:
		drop.append([num,column])
		counter_q+=1
	else:
		one = False
		for i in df_drop[column].unique():
			if len(i) != 1:
				one=True
		if one == True:
			drop.append([num,column])
			counter_len +=1

print('dropped because of length:', counter_len)
print('dropped because of being the same as ref:',counter_uniq)
print('dropped because of to low dp:',counter_dp)
print('dropped because of to low QUAL',counter_q)

#remove snps from dataframe that did not pass the filter
retract = 0
for i in drop:	
	df = df.drop(i[1], axis=1)
	del parchement_snp[i[0]-retract]
	del parchement_dp[i[0]-retract]
	del parchement_qd[i[0]-retract]
	del parchement_qual[i[0]-retract]
	retract +=1

print(df[0:5])


#create biosample column
bs_list = list(bs[1])
bs_list.insert(0, 'Reference cow')
df['biosample'] = bs_list
#create breed column
br_list = list(bs[2])
br_list.insert(0, 'Reference cow')
df['breed'] = br_list

#need to do something to get the name there
parchement_snp.append('Parchment  ')
parchement_snp.append('Parchment  ')
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

print('transform to nexus format')
list_of_lists = []
indexs = []
rows = []

for i in range(len(df)):
	seq = ''.join(list(df.iloc[i])[:-2])
	indexs.append(i)
	rows.append(seq)

#nexus file format can only have maximum of 60 characters on a line
#therefore incase the length is over 60 split it in different parts
frac, whole = math.modf(len(rows[0])/60)
last_bit = []

if int(whole) > 0:
	for i in range(int(whole)):
		list_60 = []

		for seq in rows:
			list_60.append(seq[(i * 60):((i+1) * 60)])

			if i + 1 == int(whole):
				last_bit.append(seq[((i+1) * 60):len(seq)])

		list_of_lists.append(list_60)
		if len(last_bit) != 0:
			list_of_lists.append(last_bit)

else:
	list_60 = []
	for seq in rows:
		list_60.append(seq)
	list_of_lists.append(list_60)

print('write to nexus file')
with open(args.nexus, 'w') as f:
	f.write('#nexus\n\n')
	f.write('BEGIN Characters;\n')
	f.write('DIMENSIONS ntax=%i nchar=%i;\n' % (len(list_of_lists[0]), len(rows[0])))
	f.write('FORMAT\n')
	f.write("\tdatatype='dna' missing=? gap=- symbols=\"atgc\" labels=left transpose=no interleave=yes;\n")
	f.write('MATRIX\n') 
	for block in list_of_lists:
		for idx in range(len(indexs)):
			#put indexs[idx] in first place for biosample names
			f.write( '\'%s\'\t%s\n' % (idx, block[idx].lower()))
		f.write('\n')
	f.write(';\nEND; [Characters]')
