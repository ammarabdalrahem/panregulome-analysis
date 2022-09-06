#! /usr/bin/bash

#BUSCO v5.2.2

    # Perform BUSCO analysis on each assembly file
    ls -1 out/uncom_data/*.fa| 
    sort -u |
    while read i
    do
        NAME_base=$(basename $i)
        echo "$NAME_base"
        busco -o ${NAME_base} -i $i -l poales -m genome --out_path out/busco -f 
    done
