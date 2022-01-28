#! /usr/bin/bash

#konw your g++ version
#g++ -v (for me is 9) and modify it in makefile

# creat masking_repeats directory
mkdir -p masking_repeats

#move to masking_repeats directory
cd masking_repeats

# download Scripting analyses of genomes in Ensembl Plants by cloning
git clone https://github.com/Ensembl/plant-scripts.git

#move plant-scripts directory, that creat after cloning 
cd plant-scripts

#installation steps for Repeat masking and annotation tools
make install_repeats #need a local machine password
make install
make install install_repeats
cd repeats


# run repeat masking and annotation tools for all fasta files
ls -1 ../../../*.fa| # move to work directory that contain fasta format 
    sort -u | #sort all files by alphabetical arrangement
    while read i
    do
        NAME_base=$(basename $i)
        echo "$NAME_base"
        ./Red2Ensembl.py ../../../${NAME_base}  ${NAME_base}_file \
--msk_file ${NAME_base}.sm.fna --bed_file ${NAME_base}.bed --cor 4 && ./AnnotRedRepeats.py \
../files/nrTEplantsJune2020.fna ${NAME_base}_file --bed_file ${NAME_base}.nrTEplants.bed --cor 4

    done


