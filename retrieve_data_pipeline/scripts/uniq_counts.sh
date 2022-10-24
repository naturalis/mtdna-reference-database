#!/bin/bash
#Fourth step (optional)
#Takes all the output from seperate_lines.py, sorts and gets the counts of the uniq breeds
#in the files. Only usefull to see how many you have of each breed.
cat result/breeds/*.txt | sort | uniq -c > result/all_breed_counts.txt
