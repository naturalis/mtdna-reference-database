# Scripts

### extract_bio_sample.py

The extract_bio_sample.py script takes a total of two commandline arguments:

    -i  --filename 
                      Path to either a SraRunTable or ENA file report.
    -o  --outfile 
                      The output file to write the biosamples to

The script only extracts the biosample IDs from the file, saves and passes these then on to the next script.

### seperate_lines.py

The seperate_lines.py script takes a total of two arguments:

    -i  --filename 
                      Path to input file containing descriptions for the biosample IDs.
    -o  --outfile 
                      The output file to write the breed of each biosample to.

The script takes a file containing the descriptions retrieved from NCBI as it's input. Based on which common starting identifier the biosamples have X amount of lines are retrieved from the file. From these lines the name of the breed is extracted using a regular expression.

### sample_breed.py

The sample_breed.py script takes the following commandline arguments:

    -i  --filename 
                      Path to the original input files for the pipeline.
    -o  --outfile 
                      The output file to write the breed of each biosample to.
    -d --descriptions
    				          The description files genereated by the seperate_lines.py script.
    -b --breedfile
    				          The path to the 'all_breed_counts.txt' file.

The script takes the breedfile, descriptions and the original input files and reformates all the usefull information into one CSV files.