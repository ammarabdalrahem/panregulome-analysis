#!/usr/bin/env python3

# Import Library
import statistics
import sys
import glob

#how we extract mash distance from mash output
#for i in $(ls *_mash) ; do cat $i| perl -pe 's/\t/,/g' |cut -d "," -f 3 |paste -sd , > $i.dist ; done  #to creat file

input_file=sys.argv[1] #input file name
number_list=[]


for filename in glob.glob(input_file, recursive=True):
    file_name=str(filename).split("/")[-1].split("_mash")[0] # file name
    with open(filename) as file:
        for line in file.read().split(","): # read file as lines split it by ,
            number_list.append(float(line)) # collect data as float data in list 
result =statistics.median(number_list)      # calculate median for number list
print (file_name,",","Med_mash_dist,",result)

#how we calculate the median script for all cluster and save it 
#for i in $(ls *.dist) ; do ./../../../src/scripts/median.py $i ; done > ../../../tables/total_med_mash_dist_pr.csv
#for i in $(ls *.dist) ; do ./../../../src/scripts/median.py $i ; done > ../../../tables/total_med_mash_dist_cds.csv
