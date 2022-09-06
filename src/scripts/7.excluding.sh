#! /usr/bin/bash

#save the R output in txt file
Rscript src/scripts/boxplot.R > r_output.txt

#extract the two final line with ecotypes name | cut it | replace "" with * as regx expression for loop
cat r_output.txt|tail -n 2|cut -d " " -f 3-15|sed 's/"/*/g' > excluded.txt

#note: outlier more than two line modify it to be universal??
#take a list from outlier of boxplot in R to sure good extraction
echo "$(<excluded.txt)"

#clear data
#move excluded_ecotypes to out directory in uncompress ?? how can i do it
for file in $(echo "$(<excluded.txt)"); do mv "$file" ./out/excluded done
mv excluded_list.txt r_output.txt out/txt



