#! /usr/bin/bash

#install csvjoin by python
pip install csvkit

output_ditectory=$1 #work_directory 

# creat result directory
mkdir -p result
 
#copy  all short summary result of assemblies to result directory 
cp ./*_busco/short_summary.* ./result/ 
#moving to result directory 
cd  "result"

# extract the busco score of completeness and create a column
ls -1 short_summary*.txt | 
    while read i
    do
    cat $i |sed -n '3p'| rev |cut -d "/" -f 1|rev|tr -s " " && cat $i|
    perl -ne 'if(/C:([^\[]+)/){ print $1; }'|
    awk '{print $0","}' 
    done |paste -sd "\t\n"|sort -u| awk '{print $2}'|sed 's/.$//' |sed 's/.$//'> busco_table.csv
mv busco_table.csv ..

# moving to work directory
cd ..

    
#create a table for busco column and stats result to create an assembly assessment summary table

paste -d , assessment_table_.csv busco_table.csv > final_table.csv
#creat headers
sed -i -e '1i Name,total_length,N_count,Gaps,N_genes,BUSCO_complete %' final_table.csv

#clear row names 
cat final_table.csv | awk '{gsub(/Bdistachyon/,"",$1)}1' | awk '{gsub(/_v1.fa/,"",$1)}1' > quality_table.csv

