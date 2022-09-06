#! /usr/bin/bash

#go to bed files directory
cd masking_repeats/plant-scripts/repeats/

# cut TE families after # from nrTEplants.bed files and add file -ecotype- name to columns 
ls -1 *.nrTEplants.bed | 
    sort -u |
    while read i
    do
        perl -lane '$f=(split(/#/,$F[3]))[1];print "$f\t$ARGV"' $i  
            done | tr "\\t" ", " > te_ann_0.csv   # make comma as separator 

#rearrange , put ecotype name as first column
cat te_ann_0.csv|awk -F, '{ print $2","$1 }'  > te_ann_1.csv

#filter name of ecotype to be easy to read |and replace /  by comma as separator| remove third col (superfamily)
cat te_ann_1.csv| awk '{gsub(/Bdistachyon/,"",$1)}1' | awk '{gsub(/_v1.fa.nrTEplants.bed/,"",$1)}1'|tr '/' ',n'|cut -d, -f3 --complement> te_ann_2.csv 

#measure the length(bp) for each TE order in nrTEplants.bed files save it as column 
ls -1 *.nrTEplants.bed | 
    sort -u |
    while read i
    do
        cat $i |awk '{ print $0, $3 - $2 }'| cut -d " " -f 2
            done | tr "\\t" ","  > te_ann_3.csv


#merge two tables for final table 
paste -d , te_ann_2.csv te_ann_3.csv > te_ann.csv 
#cat te_ann_4.csv|sed 's/^,/NULL,/; :a;s/,,/,NA,/g;ta' > te_ann.csv

#add header to table 
sed -i -e '1i Ecotype,Order,Length' te_ann.csv

#reomve files 

rm te_ann_0.csv
rm te_ann_1.csv
rm te_ann_2.csv
rm te_ann_3.csv

#move final file to work directory 
mv te_ann.csv ./../../..
