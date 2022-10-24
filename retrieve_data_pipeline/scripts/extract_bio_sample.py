#First step
#Extract all the BioSample ids from the file so we can do some bash code to see counts and the uniq breeds.
#output files are biosamples/biosamples_SRA and biosamples/biosamplesEVA.

import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='Retrieve BioSample information for SRA or EVA data from the NCBI server')
parser.add_argument('-i', dest="filename", required=True, help='Input file being either a SraRunTable (SRA) or a filereport.tsv (EVA) containing the BioSample accesion IDs.')
parser.add_argument('-o', dest="outfile", required=True, help='Name for the output file.')

arg = parser.parse_args()

#Run table files from SRA always start with "Sra" and files from the EVA start with 'file'
#files from EVA need their name changed to "SraRunTablefile[something]" otherwise the pipeline
#won't process the file
if "fil" in arg.filename:
	data = pd.read_csv(arg.filename, sep="\t", header=0)

	column = data["sample_accession"]

elif "Sra" in arg.filename :
	data = pd.read_csv(arg.filename, sep=",", header=0)

	column = data["BioSample"]
#Warning incase its not a file from SRA or EVA and just a reminder if the file just contains run and biosample IDs
else:
	print("File doens't seem to be a SraRunTable or filereport from SRA or EVA respectively. It is therefore assumed the input file contains just run and BioSample ID's. If this is not the case the process will not work and crash somewhere.")

with open(arg.outfile, 'a') as f:
	try:
		for i in column:
			f.write("%s\n" % i)
	except NameError:
		with open(arg.filename, 'r') as f2:
			for i in f2:
				f.write("%s\n" % i)