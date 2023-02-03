#Fourth step, takes the output from seperate_lines.py, the original input and a 
#file containing already known breeds
#BioSample IDs and run IDs

# -i, The SRA and EVA files that were used as the original input or the ID file
# -d, The files containing the retrieved descriptions from the NCBI server.
# -b, The breed file that has had the manual country and direction added to it, otherwise use the backup
# -o, The output file name.

import itertools
import pandas as pd 
import regex as re 
import argparse
import json

parser = argparse.ArgumentParser(description='Reformats all the needed info in one output file.')
parser.add_argument('-i', nargs='+', dest="filenames", required=True, help='The input files.')
parser.add_argument('-d', nargs='+', dest="descriptions", required=True, help='The description files.')
parser.add_argument('-b', dest="breedfile", required=True, help='all_breed_counts.txt. Contains the unique breeds in the ')
parser.add_argument('-o', dest="outfile", required=True, help='Name for the output file.')

arg = parser.parse_args()

count_samples = 0
directions = []
#arg input -b "result/backup_dont_touch/backup_species_country.txt"
df_breeds = pd.read_csv(arg.breedfile, sep = "    ", names=['breed','country'])

for i in range(len(df_breeds['breed'])):
	number = re.match("[0-9]{1,3}", df_breeds.iloc[i]['breed'])
	df_breeds.iloc[i]["breed"] = df_breeds.iloc[i]['breed'][number.end()+1:]
	#split the country and direction 
	if len(df_breeds.iloc[i]['country']) != 2:
		try:
			split_country = df_breeds.iloc[i]['country'].split('-')
			directions.append(split_country[1])
			df_breeds.iloc[i]['country'] = split_country[0]
		except IndexError:
			if len(df_breeds.iloc[i]['country']) != 5:
				directions.append("Unspecified: Doesn't conform to standard")
			else:
				directions.append("mixed")
	else:
		directions.append("WHOLE")	

df_breeds['direction'] = directions
print(df_breeds[0:5])

parts = []

#this gets all the run and biosample ids in the input files
#These are the original input files (SraRunTable, file) -i input/SraRunTable.txt input/filereport_analysis_PRJEB42783_tsv.txt
for file in arg.filenames:
	if 'fil' in file:
		data = pd.read_csv(file, sep="\t", header=0).rename(columns={"sample_accession": 'BioSample',"submitted_ftp": 'Run'})
		for i in range(len(data['Run'])):
			split_link = data['Run'].iloc[i].split(';')
			data['Run'].iloc[i] = '%s,%s' % (split_link[34],split_link[54])
	else:
		data = pd.read_csv(file, sep=",", header=0)

	columns = data[["BioSample", "Run"]]
	parts.append(columns)

full_df = pd.concat(parts, ignore_index = True)
print(full_df)
des_parts = []

#these are the description files -d result/biosample_descriptions/<code>_descriptions.txt
for file in arg.descriptions:
	#put file in df
	df_des = pd.read_csv(file, sep = '\t', names = ['breed', 'biosample', 'experiment_id'])
	#put df in list
	des_parts.append(df_des)

#concat df
df_des = pd.concat(des_parts, ignore_index=True)
print(df_des)

#when you rerun everything from the start you might run into errors because of manual changes made to the all_breed_counts.txt, therefore the backup exists.  

new_columns = []
counter = 0

with open("result/issues.txt", 'w') as f:
	for i in range(len(df_des['breed'])):
		columns = []
		try:
			columns.append(df_breeds["country"].loc[df_breeds["breed"] == df_des['breed'].iloc[i]].iloc[0])
			columns.append(df_breeds["direction"].loc[df_breeds["breed"] == df_des['breed'].iloc[i]].iloc[0])
			columns.insert(0,full_df["Run"].loc[full_df['BioSample'] == df_des['biosample'].iloc[i]].iloc[0])
			new_columns.append(columns)
		except IndexError:
			counter+=1
			f.write("%s, %s\n" % (df_des['breed'].iloc[i], df_des['biosample'].iloc[i]))
			new_columns.append(['?', '?', '?'])

print("\n%s samples were not processed correctly, see \"issues.txt\" for the breed and biosample id.\n" % counter)

df_des[['run', 'country', 'direction']] = new_columns

#create a dictionary containing the counts of which country and which direction the breeds are from.
country_dict = {}

for i in range(len(df_des)):
	code = df_des["country"].iloc[i]
	direction = df_des['direction'].iloc[i]
	if code in country_dict:
		if direction in country_dict[code]:
			country_dict[code][direction] += 1
		else:
			country_dict[code][direction] = 1
	else:
		country_dict[code] = {}
		country_dict[code][direction] = 1

for key in country_dict:
	print(key)
	for key2 in country_dict[key]:
		print("\t" + key2 + ": ", country_dict[key][key2])

with open('result/breed_spread.txt', 'w') as f:
     f.write(json.dumps(country_dict))
     
#arg.output -o result/cow_info.csv
df_des.to_csv(arg.outfile, sep='\t', index=False)


#country codes in the backup, current count in folder uses ISO Alpha-2 codes
# SC = SCotland
# SW = SWiss
# IT = ITaly
# FR = FRance
# DE = DEnmark
# DU = Duitsland
# AU = AUstria
# SL = SLovenie
# BE = BElgium
# UK = United Kingdom
# AM = AMerica
# AS = AuStralia
# EU = EUrope
# AF = AFrica
# JP = JaPan
# IN = INdia
# ?  = NAN
# TU = TUrkey
# KO = KOsovo
# MO = MOngolia
# CH = CHina
# SL = SLovania
# PA = PAragua
# UR = URugua
# GR = GReece
# AL = ALbania
# RU = RUsia
# PO = POrtugal
# SE = SErbia
# ET = EThiopia
# EU = EUrope
# MA = MAlta
# TI = TIbet
# PA = PAkistan
# CO = COlombia
# PL = PoLand
# KR = KoRea
# IR = IReland
# IQ = IraQ
# IA = IrAn
# HU = HUngaria
# FN = FiNland
# AN = ANcestor
# UA = UkrAine
# BA = Baltic
# EG = EGypt
# AZ = AZerbeijan
# AR = ARmania
# GE = GEorgia
# CR = CRoatia
# BU = BUlgaria
# MC = MaCadonia

#Directions:
#C = Center
#WHOLE = whole country
#other cardinal directions

#European list:
#[SC, SW, IT, FR, DE, DU, SL, BE, UK, EU, SL, GR, AL, PO, SE, MA, PL, IR, HU, FN, AN, UA, AZ, AR, GE, CR, BU, MC]

