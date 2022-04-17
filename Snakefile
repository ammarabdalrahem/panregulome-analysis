## Snakemake - panregulome-analysis
##
## @AbdalrahemAmmar

rule all:
    input:
        obtain_data_out = "out/data/",
        uncompressed_data_out = "out/uncom_data/",
        asses_table = "tables/assessment_table.csv",
        BUSCO_run = "out/busco/",
        quality_table = "tables/quality_table.csv",
        q_boxplot = "figures/quality_box.eps",
        excluded = ("out/txt/poor_quality.txt"),
        repeats_data = "out/masking_repeats",
        repeats_table = "tables/repeats_analysis.csv"


        


# --- Obtain genome sequences --- #
username = "login=mohamedabdelfadeel1994@gmail.com"
password = "password=Rido@1994"


rule obtain_data:
  input:
    script = "src/scripts/data_retrieve.sh",
    ecotypes = "src/seq.new.txt"
  
  output:
    data = directory("out/data/")
  shell:
    """
    cd {output.data}
    while read file; do bash ../../{input.script} {username} {password} $file; done < ../../{input.ecotypes}
    """
    
    
rule uncompressed_files:
  input:
    data = "out/data/"
  output:
    uncompressed_data = directory("out/uncom_data/")
  shell:
    """
    for file in $(ls {input.data}/*gz | cut -c 10- | rev | cut -c 4- | rev); do gunzip -k {input.data}/$file.gz; cat {input.data}/$file > {output.uncompressed_data}/$file; done &&
    cd out/uncom_data/
    mv BdistachyonBd21v2_1_283_v2.0.fa BdistachyonBd21v2_283_v2.fa
    mv BdistachyonBd21v2_1_283_Bd21v2.1.gene.gff3 BdistachyonBd21v2_283_v2.Bd21.1.gene.gff3
    mv BdistachyonBd21v2_1_283_Bd21v2.1.cds_primaryTranscriptOnly.fa BdistachyonBd21v2_283_v2.Bd21.1.cds_primaryTranscriptOnly.fa
    for f in *Only.fa ; do fnew=`echo $f|cut --complement -d '.' -f 2,3`; mv $f $fnew ; done 
    for f in *.gff3 ; do fnew=`echo $f|cut --complement -d '.' -f 2-4`; mv $f $fnew ; done
    """

# --- Quality control of the assemblies --- #

rule assembly_assessment:
  input:
    data= "out/uncom_data/"
    # install assembly stats

  output:
    asses_table = "tables/assessment_table.csv"
 
  shell:
    """
    #analysis stats for assembly files
    #assembly-stats get ssembly statistics from FASTA and FASTQ files.
    assembly-stats -u {input.data}/*_{{v1.fa,v2.fa}}|cut -f 1-2,7,8> stats.csv

    #convert output to CSV file
    cat stats.csv |tr "\\t" "," > stats_table.csv

    #loop for count and cut genes numbers
    for i in {input.data}/*.gff3; do
    sed -n '3p' $i | cut -d ' ' -f 3,4 && cat $i| grep -v "#" | cut -f3|sort | uniq -c|grep "gene" |grep  -o '[[:digit:]]*'
    done | paste -sd "\t\n" > genes.csv

    #cut only gene numbers
    cut -d" " -f2 genes.csv |awk '{{print $2}}' > genes_table.csv

    #join two tables
    paste -d , stats_table.csv genes_table.csv|perl -l -F/ -e 'print "$F[3]"'  > {output.asses_table}

    #clear data
    rm genes_table.csv genes.csv stats_table.csv stats.csv
    """


rule BUSCO_run:
  input:
    data= "out/uncom_data"
    # install BUSCO

  output:
    BUSCO_run = directory("out/busco/")

  shell:
    """
    #BUSCO v5.2.2

    # Perform BUSCO analysis on each assembly file
    ls -1 {input.data}/*_{{v1.fa,v2.fa}}| 
    sort -u |
    while read i
    do
        NAME_base=$(basename $i)
        echo "$NAME_base"
        busco -o ${{NAME_base}} -i $i -l poales -m genome --out_path {output.BUSCO_run} -f 
    done && rm -r busco_*
    """

rule quality_table:
  input:
    data = "out/busco/",
    asses_table = "tables/assessment_table.csv"

  output:
    quality_table = "tables/quality_table.csv"

  shell:
    """
    # extract the busco score of completeness and create a column
    ls -1  {input.data}*_{{v1.fa,v2.fa}}/short_summary.*|sort -u| 
    while read i
    do
    cat $i |sed -n '3p'| rev |cut -d "/" -f 1|rev|tr -s " " && cat $i|
    perl -ne 'if(/C:([^\[]+)/){{ print $1; }}'|
    awk '{{print $0","}}' 
    done |paste -sd "\t\n"|sort -u| awk '{{print $2}}'|sed 's/.$//' |sed 's/.$//'> tables/busco_table.csv

    #create a table for busco column and stats result 
    paste -d , {input.asses_table} tables/busco_table.csv > tables/final_table.csv

    #creat headers
    sed -i -e '1i Name,total_length,N_count,Gaps,N_genes,BUSCO_complete %' tables/final_table.csv

    #clear row names 
    cat tables/final_table.csv | awk '{{gsub(/Bdistachyon/,"",$1)}}1' | awk '{{gsub(/_v1.fa/,"",$1)}}1'| 
    awk '{{gsub(/_v2.fa/,"",$1)}}1'> {output.quality_table}

    rm tables/busco_table.csv tables/final_table.csv
    
    """

rule assemblies_evaluation:
  input:
    script = "src/scripts/boxplot.R",
    quality_table = "tables/quality_table.csv"

  output:
    q_boxplot = "figures/quality_box.eps"
  
  shell:
    "Rscript {input.script}"


rule excluded_poor_quality:
  input:
    script = "src/scripts/boxplot.R"

  output:
    excluded = ("out/txt/poor_quality.txt")

  
  shell:
    """    
    #save the R output in txt file
    Rscript {input.script} > r_output.txt

    #extract the two final line with ecotypes name | cut it | replace "" with * as regx expression for loop
    cat r_output.txt|tail -n 2|cut -d " " -f 3-15|sed 's/"/*/g' > poor_quality.txt

    #note: outlier more than two line modify it to be universal??
    #take a list from outlier of boxplot in R to sure good extraction
    echo "$(<excluded.txt)"

    #clear data
    #move excluded_ecotypes  
    for file in $(echo "$(<excluded.txt)"); do mv "$file" {output.excluded}; done
    mv poor_quality.txt {output.excluded}
    rm r_output.txt 
    """

    
# --- Repeat masking --- #

rule repeat_masking_annotation:
  input:
    script = "src/scripts/msk_ann.sh"

  output:
    repeats_data = directory ("out/masking_repeats")
  
  shell:
    """
    #konw your g++ version
    #g++ -v (for me is 9) and modify it in makefile

    # creat masking_repeats directory
    mkdir -p out/masking_repeats

    #move to masking_repeats directory
    cd out/masking_repeats

    # download Scripting analyses of genomes in Ensembl Plants by cloning
    git clone https://github.com/Ensembl/plant-scripts.git

    #move plant-scripts directory, that creat after cloning 
    cd plant-scripts
    #installation steps for Repeat masking and annotation tools
    make install_repeats #need a local machine password how m
    make install
    cd repeats
    #REPEATS=out/plant-scripts/repeats


    # run repeat masking and annotation tools for all fasta files
    ls -1 ../../../uncom_data/*_{{v2.fa,v1.fa}}| # move to work directory that contain fasta format 
    sort -u | #sort all files by alphabetical arrangement
    while read i
    do
        NAME_base=$(basename $i)
        echo "$NAME_base"
        ./Red2Ensembl.py ../../../uncom_data/${{NAME_base}}  ${{NAME_base}}_file \
      --msk_file ${{NAME_base}}.sm.fna --bed_file ${{NAME_base}}.bed --cor 4 && ./AnnotRedRepeats.py \
      ../files/nrTEplantsJune2020.fna ${{NAME_base}}_file --bed_file ${{NAME_base}}.nrTEplants.bed --cor 4

    done

    mv *.bed *_file ../../../repeats
    """

rule repeats_length:
  input:
    data = "out/masking_repeats/plant-scripts/repeats/"

  output:
    repeats_table = "tables/repeats_analysis.csv"
  
  shell:
    """
    #count total length of repeats for each ecotypes
    ls -1 {input.data}*.fa.bed | 
        sort -u |
        while read i
        do
            cat $i |awk '{SUM += $3-$2} END {print SUM}' # by subtracting two columns of coordinates for length the summation all for total length
                done | tr "\\t" ","  > repeats_length.csv

    cat repeats_length.csv

    #count the total length of genome for each ecotypes
    assembly-stats -u Bdistachyon*.fa|cut -f 1-2|tr "\\t" ","|awk '{gsub(/Bdistachyon/,"",$1)}1' | awk '{gsub(/_v1.fa/,"",$1)}1' > stats_40.csv

    # join genome total length and total length of repeats in CSV file
    paste -d , stats_40.csv repeats_length.csv > tables/repeats_analysis.csv 

    #add head 
    sed -i '1i Ecotypes,Genome.size,Repeated' tables/repeats_analysis.csv

    #clean data
    rm stats_40.csv repeats_length.csv
    """ 



rule repeats_visualization:
  input:
    data = "out/masking_repeats/plant-scripts/repeats/"

#  output:
#    data = "src/work/*.csv"
  
  shell:
    "Rscript {input.script}"

'''
rule annotation_analysis:
  input:
    script = "src/scripts/te_ann.sh"

#  output:
#    data = "src/work/*.csv"
  
  shell:
    "bash {input.script}"

rule annotation_visualization:
  input:
    script = "src/scripts/te_orders.R"

#  output:
#    data = "src/work/*.csv"
  
  shell:
    "Rscript {input.script}"


# --- Extraction of proximal promoter sequences --- #

rule promoters_extraction:
  input:
    script = "src/scripts/promoters_extraction.sh"

#  output:
#    data = "src/work/*.csv"
  
  shell:
    "bash {input.script}"





#snakemake --rulegraph| dot -Tpng > rulegraph.png
#snakemake --dag | dot -Tpng > ruledag.png

import os, glob 

samples = []

for file in glob.glob("src/data/*fa"):
    name = file.split("/")[-1].split(".")[0]
    #print(name)
    if name not in samples:
        samples.append(name)
        #print(samples)
   

'''