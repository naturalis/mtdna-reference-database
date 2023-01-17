# mtdna-reference-database

Reference database of panels of mitochondrial genomes for population-level aDNA identification

Table of file codes related to pieces of parchment:

|Code|Parchment|
|----|---------|
|P1_S7|Object 1|
|P2_S8|Object 2|
|P3_S9|Object 3|
|P5_S6|Object 5|
|P195_S1|195|
|P700_S2|700|
|P786_S5|786|
|P78404_S3|784-04|
|P78406_S4|784-06|

This repository consits of the code, pipelines, databases and documentation for recreating the results of the parchment project. These results consist of the reference Bos Taurus database and the vcf files for all the samples of parchment.The following
subsections are here:

- [Documentation](doc), about packages and tools needed for everything in the repository.
- [Cow_metadata_pipeline](cow_metadata_pipeline), containing the files for the pipeline, including the scripts, test data and documentation 
- [Database](database), containing the scripts needed for transforming the data to the correct format, database schema, testdata and documentation.
- [Variant_calling_pipeline](variant_calling_pipeline), containing the files for the variant calling pipeline, including the scripts, test data and documentation 

reserved zenodo doi:

10.5281/zenodo.7371345

