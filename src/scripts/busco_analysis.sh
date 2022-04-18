#! /usr/bin/bash

install csvjoin by python
pip install csvkit
 

# extract the busco score of completeness and create a column
ls -1  out/busco/*.fa/short_summary.*|
    while read i
    do
    cat $i |sed -n '3p'| rev |cut -d "/" -f 1|rev|tr -s " " && cat $i|
    perl -ne 'if(/C:([^\[]+)/){ print $1; }'|
    awk '{print $0","}' 
    done |paste -sd "\t\n"|sort -u| awk '{print $2}'|sed 's/.$//' |sed 's/.$//'> tables/busco_table.csv

#create a table for busco column and stats result to create an assembly assessment summary table

paste -d , tables/assessment_table.csv busco_table.csv > final_table.csv
#creat headers
sed -i -e '1i Name,total_length,N_count,Gaps,N_genes,BUSCO_complete %' final_table.csv

#clear row names 
cat final_table.csv | awk '{gsub(/Bdistachyon/,"",$1)}1' | awk '{gsub(/_v1.fa/,"",$1)}1' > quality_table.csv

# clear data
#rm busco_table.csv assessment_table.csv 
