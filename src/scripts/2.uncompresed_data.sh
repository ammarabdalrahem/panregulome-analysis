#!usr/bin/env bash

rm BrachyPan.file_list cookies

 #change the directory
gunzip *.gz
mkdir annotation
mkdir assembly
#find . -name '*.gff3' -exec mv {} ./annotation/ ';'
#find . -name '*.fa' -exec mv {} ./assembly/ ';'

# to run  bash /home/ammar/ammar/bash.script/uncompresed_data.sh 
