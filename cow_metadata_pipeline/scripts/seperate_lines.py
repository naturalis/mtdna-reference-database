#Third step, this script takes the output from get_description.py
#loop through file and retrieve only breed name to store in a file.

#input can be seen in if statements
#output will be result/breeds/[code]_breeds

import itertools
import pandas as pd
import argparse
import regex as re

parser = argparse.ArgumentParser(description='Extract the breed of the cow from the retrieved BioSample information.')
parser.add_argument('-i', dest="filename", required=True, help='Input file being the file containing the feteched BioSample descriptions.')
parser.add_argument('-o', dest="outfile", required=True, help='Name for the output file.')

arg = parser.parse_args()
b_pattern = re.compile('\D{3,5}\d{7,8}')
r_pattern = re.compile('\D{3}\d{6,8}')

if 'None' in arg.filename:
	exit()
#opens the input file and only retrieves x amount of lines. The correct line that should 
#contain the breed name is then searched over with a regular expression to retrieve the full
#breed name and then stored in a file.
with open(arg.filename, 'r') as f:
	with open(arg.outfile, "w") as f2:
		if "SAMEA" in arg.filename:
			for l1,l2,l3,l4,l5,l6 in itertools.zip_longest(*[f]*6):
				breed = re.search("[A-Za-zè -.]*$", l1).group()[1:]
				if len(breed) > 40:
					continue
				if breed == "e":
					breed = 'Evolène'
				if breed == "":
					if bool(re.search(r'\d', l4)) or l4.strip("\n")[-1] == "\"":
						continue
					else:
						breed = l4.strip("\n")
				else:
					pass

				accessions = l2.strip("\n").split(";")
				biosample = b_pattern.search(accessions[0]).group()
				run = r_pattern.search(accessions[1]).group()

				f2.write("%s\t%s\t%s\n" % (breed, biosample, run))

		elif "SAMN" in arg.filename or "SAMD" in arg.filename:
			for line1,line2 in itertools.zip_longest(*[f]*2):
				if line2[:5] != "    /":
					try:
						breed = re.search('\".*\"', line1).group()
						accessions = line2.strip("\n").split(";")
						biosample = b_pattern.search(accessions[0]).group().strip(" ")
						run = r_pattern.search(accessions[2]).group()
					except:
						continue
				else:
					try:
						breed = re.search('\".*\"', line2).group()
						accessions = line1.strip("\n").split(";")
						biosample = b_pattern.search(accessions[0]).group().strip(" ")
						run = r_pattern.search(accessions[-1]).group()
					except:
						continue
				
				f2.write("%s\t%s\t%s\n" % (breed.strip("\""), biosample, run))


		else:
			print("Will try to retrieve the breed, but will likely crash as the place of the breed name is not know for this file: %s" % arg.filename)
			#loop over every line.
			#for line in f:
			#see if the word "breed" is in the line

			#if so then do a regex and save
			#do this for run ID(if applicable) aswell
			#else keep going to the next line
			exit()
