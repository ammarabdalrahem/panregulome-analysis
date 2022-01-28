#! /usr/bin/bash

#BUSCO v5.2.2

#define the work directory
ASSEMBLY_FOLDER=$1

#echo "define the directory"
# numer of cores
NUM_CORES=9

# Perform BUSCO analysis on each assembly file
ls -1 ${ASSEMBLY_FOLDER}/*.fa| 
    sort -u |
    while read i
    do
        NAME_base=$(basename $i)
        echo "$NAME_base"
        busco -o ${NAME_base}_busco -i $i -l poales -m genome -f -c $NUM_CORES 
    done

