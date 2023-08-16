import pandas as pd
import sys

# Ensure we have a command-line argument
if len(sys.argv) < 2:
    sys.stderr.write("Please provide the path to the TSV file as an argument.\n")
    sys.exit(1)

# Read the TSV file
data = pd.read_csv(sys.argv[1], sep='\t')

# List of unique Superregions
superregions = data['Superregion'].unique()

# Initialize a new dataframe
output = pd.DataFrame()
output['SampleID'] = data['SampleID']

# Add columns for each Superregion
for superregion in superregions:
    output[superregion] = (data['Superregion'] == superregion).astype(int)

# Write the result to stdout
output.to_csv(sys.stdout, sep='\t', index=False)
