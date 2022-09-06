#! /usr/bin/bash

#count total length of repeats for each ecotypes
ls -1 ./masking_repeats/plant-scripts/repeats/*.fa.bed | 
    sort -u |
    while read i
    do
        cat $i |awk '{SUM += $3-$2} END {print SUM}' # by subtracting two columns of coordinates for length the summation all for total length
            done | tr "\\t" ","  > repeats_length.csv

cat repeats_length.csv

#count the total length of genome for each ecotypes
assembly-stats -u Bdistachyon*.fa|cut -f 1-2|tr "\\t" ","|awk '{gsub(/Bdistachyon/,"",$1)}1' | awk '{gsub(/_v1.fa/,"",$1)}1' > stats_40.csv

# join genome total length and total length of repeats in CSV file
paste -d , stats_40.csv repeats_length.csv > repeats_analysis.csv 

#add head 
sed -i '1i Ecotypes,Genome.size,Repeated' repeats_analysis.csv 

