#!usr/bin/bash
#assembly-stats

#analysis stats for assembly files 
#assembly-stats get ssembly statistics from FASTA and FASTQ files.
assembly-stats -u Bdistachyon*.fa|cut -f 1-2,7,8> stats.csv

#convert output to CSV file
cat stats.csv |tr "\\t" "," > stats_table.csv

#loop for count and cut genes numbers  
for i in Bdistachyon*.gff3; do
  sed -n '3p' $i | cut -d ' ' -f 3,4 && cat $i| grep -v "#" | cut -f3|sort | uniq -c|grep "gene" |grep  -o '[[:digit:]]*'
  done | paste -sd "\t\n" > genes.csv

#cut only gene numbers 
cut -d" " -f2 genes.csv |awk '{print $2}' > genes_table.csv

#join two tables  
paste -d , stats_table.csv genes_table.csv > assessment_table_.csv

#now have  5 columns next step add buscoresult for each assembly
