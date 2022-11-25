#!/bin/bash
#Second step, takes the output from extract_bio_sample.py
#this gets all the uniq starts of ids from all the biosample files
cat intermediaries/biosamples/* | sort | egrep -o "[A-Z]{3,5}" | uniq > intermediaries/unique_identifiers.txt
echo "Retrieved unique starting identifiers for the biosamples."
#when inputting multiple files get only the unique ids put them in seperate files and fetch them seperately since they have their info stored in different places
for OUTPUT in $(cat intermediaries/unique_identifiers.txt)
do
	cat intermediaries/biosamples/* | sort | uniq | grep $OUTPUT > intermediaries/seperate_biosample_ids/$OUTPUT.txt
	if [ $OUTPUT = "SAMEA" ]
	then
		efetch -input "intermediaries/seperate_biosample_ids/$OUTPUT.txt" -db biosample -format native -mode text | grep -B1 "Accession:\|Identifiers:" > result/biosample_descriptions/SAMEA_descriptions.txt
		echo "Retrieved discriptions for SAMEA biosamples"
	elif [ $OUTPUT = "SAMN" ]
	then
		efetch -input "intermediaries/seperate_biosample_ids/$OUTPUT.txt" -db biosample -format native -mode text | grep 'breed=\|Sample name=\|Identifiers:' > result/biosample_descriptions/SAMN_descriptions.txt
		echo "Retrieved discriptions for SAMN biosamples"
	elif [ $OUTPUT = "SAMD" ]
	then
		efetch -input "intermediaries/seperate_biosample_ids/$OUTPUT.txt" -db biosample -format native -mode text |  grep "sample comment=\|Identifiers:" > result/biosample_descriptions/SAMD_descriptions.txt
		echo "Retrieved discriptions for SAMD biosamples"
	else
		echo "BioSample code does not match any code for which the script was made. The descriptions of the BioSamples are therefore retrived in full, but wil likely not be processed correctly later on in the pipeline."
		efetch -input "intermediaries/seperate_biosample_ids/$OUTPUT.txt" -db biosample -format native -mode text > result/biosample_descriptions/unknown_descriptions.txt
	fi
done









