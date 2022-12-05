import sqlite3 as sq
from sqlite3 import Error
import argparse
import pandas as pd

parser = argparse.ArgumentParser(prog = 'db_stuff.py', description = 'does db stuff.')
parser.add_argument('-i', '--input', help='path to database.', required=True )
parser.add_argument('-v', '--vcf', help='path to vcf file.', required=True )  
parser.add_argument('-o', '--output', help='path to output file.', required=True )           

args = parser.parse_args()

df_vcf = pd.read_csv(args.vcf, sep='\t',  names =['#CHROM', 'POS', 'REF', 'ALT', 'QUAL', 'SAMPLE'])



try:
    connection = sq.connect(args.input)
    print("Connection to SQLite DB successful")
except Error as e:
    print(f"The error '{e}' occurred")



#get the reference table
cursor = connection.cursor()

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



col_list = []
parchement_snp = []
for col in df_all.columns:
	if col in list(df_vcf['POS']):
		col_list.append(df_all[col])
		parchement_snp.append(df_vcf.loc[df_vcf['POS'] == col]['ALT'].iloc[0])

df = pd.concat(col_list, axis=1)
df.columns = df.iloc[0]
df = df.drop(df.index[0])



# get the biosamples and save in a list, gonna be index
try:
    cursor.execute('SELECT accession_id, biosample, breed FROM accessions;')
    acc_info = cursor.fetchall()

except Error as e:
    print(f"The error '{e}' occurred")



# for i in biosample ids
bs = list(zip(*acc_info))

for i in bs[0]:
	#	get the snps
	#	fill a row in the df with the snps and where there is no snp put a '-'
	row = []
	try:
	    cursor.execute('SELECT refsite_id, alt_allele FROM snps WHERE accession_id = %i;' % (i))
	    snp_info = cursor.fetchall()
	    snp = list(zip(*snp_info))
	    for column in df.columns:
	    	if column in snp[0]:
	    		row.append(snp[1][snp[0].index(column)])
	    	else:
	    		row.append('.')
	    df.loc[len(df)] = row

	except Error as e:
	    print(f"The error '{e}' occurred")


#create biosample 
bs_list = list(bs[1])
bs_list.insert(0, 'Reference cow')			    
df['biosample'] = bs_list

#need to do something to get the name there
parchement_snp.append('Parchment P5')
df.loc[len(df)] = parchement_snp
#set biosample as index
#df = df.set_index('biosample')

#make the MSA
for i in df.columns[:-1]:
	ref_len = len(df[i][0])
	max_len = df[i].str.len().max()
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


		if len(df[i][idx]) == max_len and df[i][idx] != '*':
			idx += 1
			continue
		else:
			diffrence = max_len - len(df[i][idx])

			if df[i][idx] == '.': 

				if ref_len < max_len:
					ref_dif = max_len - ref_len
					df[i][idx] = (ref_len * '.') + (ref_dif * '-')
				else:
					df[i][idx] = df[i][idx] + (diffrence * '.')

			elif df[i][idx] == '*':
				df[i][idx] = '-' + (diffrence * '-')

			else:
				df[i][idx] = df[i][idx] + (diffrence * '-')

		idx += 1

df.to_csv(args.output, sep='\t')

