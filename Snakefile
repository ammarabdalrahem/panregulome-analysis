## Snakemake - panregulome-analysis
##
## @AbdalrahemAmmar
#from glob import glob
#lst = glob("out/cds_est_homologues/BdistachyonABR2337v1_4taxa_algOMCL_e1_/*.fna", recursive = False)
#fna_files = []
#for i in lst:
#    new_name = i.split("/")[-1]
#    fna_files.append(new_name)
#    #print(new_name)

##------------rule all ----------##

#rule all:
#    input:
#        header_extraction_out = expand(["out/uncom_data/cds/clusters_cds_ids/{file}"], file = fna_files),
#	      promoter_out = ["out/uncom_data/promoter_cluster/"],
#        coding_sequance_align = expand(["out/aln_cds/{file}"], file = fna_files),
#        promoter_sequance_align = expand(["out/aln_promoter/{file}"], file = fna_files),
#        cds_global_align_filtration =  expand(["out/aln_cds/aln_cds_filtered/{file}"], file = fna_files),
#        promoter_global_align_filtration = expand(["out/aln_promoter/aln_promoter_filtered/{file}"], file = fna_files),
#        coding_sequance_local_align = expand(["out/aln_local_cds/{file}"], file = fna_files),
#        promoter_sequance_local_align = expand(["out/aln_local_pr/{file}"], file = fna_files),
#        trimal_cds_local_align = expand(["out/aln_local_cds/aln_local_cds_trimal/{file}"], file = fna_files),
#        trimal_pr_local_align = expand(["out/aln_local_pr/aln_local_promoter_trimal/{file}"], file = fna_files),
#        trimal_cds_global_align = expand(["out/aln_cds/aln_cds_trimal/{file}"], file = fna_files),
#        trimal_pr_global_align = expand(["out/aln_promoter/aln_promoter_trimal/{file}"], file = fna_files),
#        diversity_cds_global = "out/diversity_cds_global/{file}.txt",
#        diversity_pr_global = "out/diversity_pr_global/{file}.txt",
#        diversity_cds_local = "out/diversity_cds_local/{file}.txt",
#        diversity_pr_local = "out/diversity_pr_local/{file}.txt"



# --- Create directories --- #
rule make_directories:
  shell:
    """
    mkdir -p  tables out figures
    """

# --- Obtain genome sequences --- #
username = "login=mohamedabdelfadeel1994@gmail.com"
password = "password=Rido@1994"
# Make sure your account has permission to download data

rule obtain_data:
  input:
    script = "src/scripts/data_retrieve.sh",
    ecotypes = "src/files_list.txt"
  
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
    for file in $(ls {input.data}*gz  |rev|cut -d "." -f 2- |rev|cut -d "/" -f 3-); do gunzip -k {input.data}$file.gz ; mv {input.data}$file {output.uncompressed_data} ; done
    cd out/uncom_data/
    mv BdistachyonBd21v2_1_283_v2.0.fa BdistachyonBd21v2_283_v2.fa
    mv BdistachyonBd21v2_1_283_Bd21v2.1.gene.gff3 BdistachyonBd21v2_283_v2.Bd21.1.gene.gff3
    mv BdistachyonBd21v2_1_283_Bd21v2.1.cds_primaryTranscriptOnly.fa BdistachyonBd21v2_283_v2.Bd21.1.cds_primaryTranscriptOnly.fa
    for f in *Only.fa ; do fnew=`echo $f|cut --complement -d '.' -f 2,3`; mv $f $fnew ; done 
    for f in *.gff3 ; do fnew=`echo $f|cut --complement -d '.' -f 2-4`; mv $f $fnew ; done
    mkdir genomes
    mkdir cds
    mv ./*Only.fa cds/
    mv ./*.*  genomes/
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
    mv genes_table.csv genes.csv stats_table.csv stats.csv tables/
    """

rule BUSCO_run:
  input:
    data= "out/uncom_data/"
    # install BUSCO

  output:
    BUSCO_run = directory("out/busco/"),
    BUSCO_log = directory ("out/logs/")

  shell:
    """
    #BUSCO v5.2.2
    # Perform BUSCO analysis on each assembly file
    ls -1 {input.data}*_{{v1.fa,v2.fa}}|
    sort -u |
    while read i
    do
        NAME_base=$(basename $i)
        echo "$NAME_base"
        busco -o ${{NAME_base}} -i $i -l poales -m genome --out_path {output.BUSCO_run} -c 20 -f 
    done 
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
    sed -i -e '1i Name,total_length,Ns_count,Gaps,N_genes,BUSCO_complete %' tables/final_table.csv
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
    excluded = ("out/txt/poor_quality.txt"),
    excluded_dir = directory("out/uncom_data/poor_quality/")


  shell:
    """
    #save the R output in txt file
    Rscript {input.script} > r_output.txt
    #extract the two final line with ecotypes name | cut it | replace "" with * as regx expression for loop
    cat r_output.txt|tail -n 2|cut -d " " -f 3-15|sed 's/"/*/g' > {output.excluded}
    #note: outlier more than two line modify it to be universal??
    #take a list from outlier of boxplot in R to sure good extraction
    echo "$(<{output.excluded})"
    #clear data
    #move excluded_ecotypes
    for file in $(echo $(<out/txt/poor_quality.txt) ); do mv out/uncom_data/genomes/$file {output.excluded_dir}; done
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
    mkdir -p src/scripts/masking_repeats

    #move to masking_repeats directory
    cd src/scripts/masking_repeats

    # download Scripting analyses of genomes in Ensembl Plants by cloning
    git clone https://github.com/Ensembl/plant-scripts.git

    #move plant-scripts directory, that creat after cloning
    cd plant-scripts
    #installation steps for Repeat masking and annotation tools
    make install install_repeats
    cd repeats


    # run repeat masking and annotation tools for all fasta files
    ls -1 ../../../../../out/uncom_data/genomes/*.fa| # move to work directory that contain fasta format
        sort -u | #sort all files by alphabetical arrangement
        while read i
        do
            NAME_base=$(basename $i)
            echo "$NAME_base"
            ./Red2Ensembl.py  ../../../../../out/uncom_data/genomes/${{NAME_base}}  ${{NAME_base}}_file \
            --msk_file ${{NAME_base}}.sm.fna --bed_file ${{NAME_base}}.bed --cor 16 && ./AnnotRedRepeats.py \
            ../files/nrTEplantsJune2020.fna ${{NAME_base}}_file --bed_file ${{NAME_base}}.nrTEplants.bed --cor 16

        done

    """

rule repeats_length:
  input:
    data = "src/scripts/masking_repeats/plant-scripts/repeats/"

  output:
    repeats_table = "tables/repeats_analysis.csv"
  
  shell:
    """
    #count total length of repeats for each ecotypes
    ls -1 {input.data}*.fa.bed | 
        sort -u |
        while read i
        do
            cat $i |awk '{{SUM += $3-$2}} END {{print SUM}}' # by subtracting two columns of coordinates for length the summation all for total length
                done | tr "\\t" ","  > repeats_length.csv
    cat repeats_length.csv
    #count the total length of genome for each ecotypes
    assembly-stats -u out/uncom_data/genomes/Bdistachyon*.fa|cut -f 1-2|tr "\\t" ","|awk '{{gsub(/Bdistachyon/,"",$1)}}1' | awk '{{gsub(/_v1.fa/,"",$1)}}1'|cut -d "/" -f 4 > stats_40.csv
    # join genome total length and total length of repeats in CSV file
    paste -d , stats_40.csv repeats_length.csv > tables/repeats_analysis.csv 
    #add head 
    sed -i '1i Ecotypes,Genome.size,Repeated' tables/repeats_analysis.csv
    #clean data
    rm stats_40.csv repeats_length.csv
    """ 

rule repeats_visualization:
    input:
      script = "src/scripts/repeat_plot.R"

    output:
      data = "figures/repeats_plot.eps"

    shell:
      "Rscript {input.script}"


rule annotation_analysis:
  input:
    script = "src/scripts/te_ann.sh"
  output:
    data = "tables/te_ann.csv"
  
  shell:
    """
    #go to bed files directory
    cd src/scripts/masking_repeats/plant-scripts/repeats/

    # cut TE families after # from nrTEplants.bed files and add file -ecotype- name to columns
    ls -1 *.nrTEplants.bed |
        sort -u |
        while read i
        do
            perl -lane '$f=(split(/#/,$F[3]))[1];print "$f\t$ARGV"' $i
                done | tr "\\t" ", " > te_ann_0.csv   # make comma as separator

    #rearrange , put ecotype name as first column
    cat te_ann_0.csv|awk -F, '{{ print $2","$1 }}'  > te_ann_1.csv

    #filter name of ecotype to be easy to read |and replace /  by comma as separator| remove third col (superfamily)
    cat te_ann_1.csv| awk '{{gsub(/Bdistachyon/,"",$1)}}1' | awk '{{gsub(/_v1.fa.nrTEplants.bed/,"",$1)}}1'|tr '/' ',n'|cut -d, -f3 --complement> te_ann_2.csv

    #measure the length(bp) for each TE order in nrTEplants.bed files save it as column
    ls -1 *.nrTEplants.bed |
        sort -u |
        while read i
        do
            cat $i |awk '{{ print $0, $3 - $2 }}'| cut -d " " -f 2
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
    cp te_ann.csv ./../../../../../tables/
    """

rule annotation_visualization:
  input:
    script = "src/scripts/te_orders.R"
  output:
    data = "figures/TE_orders_plot.eps"
  
  shell:
    "Rscript {input.script}"

# --- Extraction of proximal promoter sequences --- #
rule promoters_extraction:
  input:
    script = "src/scripts/promoters_extraction.sh"
  output:
    data_pr = directory ("out/uncom_data/promoter/")
  
  shell:
    "bash {input.script}"



rule get_homologues:
  #input:
    #script = "src/scripts/get_homologues/get_homologues-est.pl"
  output:
    get_homologues = "src/scripts/get_homologues/get_homologues-est.pl"

  shell:
    """
    cd src/scripts/
    git clone https://github.com/eead-csic-compbio/get_homologues.git
    cd get_homologues
    ./install.pl
    cd ../../..
    {output.get_homologues} -h
    """
    
rule gene_clustering:
  input:
    script_get_homologues = "src/scripts/get_homologues/get_homologues-est.pl"
  output:
    get_homologues_log = "out/logs/log.cds.Pfam"

  shell:
    """
    {input.script_get_homologues} -d out/uncom_data/cds -M -t 4 -e -m cluster &> out/logs/log.cds.Pfam
    cp -r cds_est_homologues/ out/
    rm -r cds_est_homologues/
    """

rule pangenome_analysis:
  input:
    script_get_homologues_cluster = "./src/scripts/get_homologues/compare_clusters.pl" #i should add -t 4 to avoide cloude
  output:
    get_homologues_cluster_log = "out/logs/log.compare_clusters.cds",
    clusters_cds = directory("out/clusters_cds/")

  shell:
    """
    {input.script_get_homologues_cluster} -d out/cds_est_homologues/BdistachyonABR2337v1_4taxa_algOMCL_e1_/ \
    -o out/clusters_cds -t 4 -m -n &> {output.get_homologues_cluster_log}
    """

rule pangenome_plot:
  input:
    script_get_homologues_plot = "./src/scripts/get_homologues/parse_pangenome_matrix.pl",
    pangenome_matrix = "out/clusters_cds/pangenome_matrix_t4.tab"
  output:
    get_homologues_plot_log = "out/logs/log.parse_pangenome_matrix.cds.t4"

  shell:
    """
    {input.script_get_homologues_plot} -m  {input.pangenome_matrix} -s -n &> {output.get_homologues_plot_log}
    """

rule header_extraction:
  input:
    data_homologues = "out/cds_est_homologues/BdistachyonABR2337v1_4taxa_algOMCL_e1_/{file}",
    script_header = "src/scripts/header_extraction.py"
  output:
    cds_ids = "out/uncom_data/cds/clusters_cds_ids/{file}",

  shell:
    """
    {input.script_header} -i {input.data_homologues} -o {output.cds_ids}
    """

rule promoter_clustering:
  input:
    cds_ids = "out/uncom_data/cds/clusters_cds_ids/",
    script_seq_extraction = "src/scripts/seq_extraction_by_id.py",
    promoter_dir = "out/uncom_data/promoter"
  output:
    promoter_cluster = directory("out/uncom_data/promoter_cluster/")
  shell:
    """
    mkdir -p {output.promoter_cluster}/
    {input.script_seq_extraction} -id {input.cds_ids} -f {input.promoter_dir} -o {output.promoter_cluster}/
    """


rule global_alignment:
  input:
    data_homologues = "out/cds_est_homologues/BdistachyonABR2337v1_4taxa_algOMCL_e1_/{file}", 
  output:
    cds_align = "out/aln_cds/{file}",
    promoter_align = "out/aln_promoter/{file}"

  shell: 
    """
    #cds
    mafft --globalpair --maxiterate 1000 {input.data_homologues} > out/aln_cds/{wildcards.file}  2> out/logs/aln_cds_log.txt  

    #promoter
    mafft --globalpair --maxiterate 1000 out/uncom_data/promoter_cluster/{wildcards.file} > out/aln_promoter/{wildcards.file} 2> out/logs/aln_pr_log.txt
    """

rule global_align_filtration:
  input:
    cds_align = "out/aln_cds/{file}",
    promoter_align = "out/aln_promoter/{file}"
  output:
    cds_align_filtration = "out/aln_cds/aln_cds_filtered/{file}",
    promoter_align_filtration = "out/aln_promoter/aln_promoter_filtered/{file}"

  params:
    filter_path = "./src/scripts/filter_global_alignment.py"
   

  shell: 
    """
    #cds
    {params.filter_path} -f {input.cds_align} -o {output.cds_align_filtration}  >> tables/global_poor_alignment_cds.csv || true
    #promoter
    {params.filter_path} -f {input.promoter_align} -o {output.promoter_align_filtration} >> tables/global_poor_alignment_pr.csv || true
    """

rule global_align_filtration_2:
  input:
    cds_align = "out/aln_cds/",
    promoter_align = "out/aln_promoter/"
  output:
    cds_align_filtration = "out/aln_cds/aln_cds_filtered2/",
    promoter_align_filtration = "out/aln_promoter/aln_promoter_filtered2/"

  params:
    filter_path = "./src/scripts/filter_global_alignment.py"
   

  shell: 
    """
    mkdir -p {output.cds_align_filtration} {output.promoter_align_filtration}

    # Parallel processing for CDS alignment
    find {input.cds_align} -maxdepth 1 -type f -name '*.fna' | parallel --jobs 10 '{params.filter_path} -f {{}} -o {output.cds_align_filtration}/{{/}} >> tables/global_poor_alignment_cds2.csv || true'

    # Parallel processing for promoter alignment
    find {input.promoter_align} -maxdepth 1 -type f -name '*.fna' | parallel --jobs 10 '{params.filter_path} -f {{}} -o {output.promoter_align_filtration}/{{/}} >> tables/global_poor_alignment_pr2.csv || true'
    """


rule local_alignment:
  input:
    data_homologues = "out/cds_est_homologues/BdistachyonABR2337v1_4taxa_algOMCL_e1_/{file}",
    promoter_cluster = "out/uncom_data/promoter_cluster/{file}"
  output:
    cds_local_align = "out/aln_local_cds/{file}",
    promoter_local_align = "out/aln_local_pr/{file}"

  shell: 
    """
    #cds
    ./src/scripts/get_homologues/annotate_cluster.pl -f {input.data_homologues} -o  {output.cds_local_align} >>  out/logs/aln_local_cds_log.txt 
    #promoter

    #promoter
    ./src/scripts/get_homologues/annotate_cluster.pl -f {input.promoter_cluster} -o {output.promoter_local_align} >>  out/logs/aln_local_pr_log.txt 

    """

rule trimal_alignment:
  input:
    cds_local_align = "out/aln_local_cds/{file}",
    promoter_local_align = "out/aln_local_pr/{file}",
    cds_align_filtration = "out/aln_cds/aln_cds_filtered/{file}",
    promoter_align_filtration = "out/aln_promoter/aln_promoter_filtered/{file}"
  output:
    trimal_cds_local_align = "out/aln_local_cds/aln_local_cds_trimal/{file}",
    trimal_pr_local_align = "out/aln_local_pr/aln_local_promoter_trimal/{file}",
    trimal_cds_global_align = "out/aln_cds/aln_cds_trimal/{file}",
    trimal_pr_global_align = "out/aln_promoter/aln_promoter_trimal/{file}"
  params:
    trimal_path = "~contrera/soft/trimal1.2/source/trimal"
  shell: 
    """
    #trimal version trimAl 1.2rev59
        
    #cds_global_filterated
    {params.trimal_path} -in {input.cds_align_filtration} -out {output.trimal_cds_global_align} \
    -automated1 >> out/logs/trimal_cds_log.txt 


    #promoter_global_filterated
    {params.trimal_path} -in {input.promoter_align_filtration} -out {output.trimal_pr_global_align} \
    -automated1 >> out/logs/trimal_pr_log.txt 

    #cds_local
    {params.trimal_path} -in {input.cds_local_align} -out {output.trimal_cds_local_align} \
    -automated1 >> out/logs/trimal_cds_local_log.txt 


    #promoter_local
    {params.trimal_path} -in {input.promoter_local_align} -out {output.trimal_pr_local_align} \
    -automated1 >> out/logs/trimal_pr_local_log.txt 

    """
  


rule nucleotide_diversity:
  input:
    trimal_cds_local_align = "out/aln_local_cds/aln_local_cds_trimal",
    trimal_pr_local_align = "out/aln_local_pr/aln_local_promoter_trimal",
    trimal_cds_global_align = "out/aln_cds/aln_cds_trimal",
    trimal_pr_global_align = "out/aln_promoter/aln_promoter_trimal"
  output:
    diversity_cds_global = "tables/diversity_cds_global.csv",
    diversity_pr_global = "tables/diversity_pr_global.csv",
    diversity_cds_local = "tables/diversity_cds_local.csv",
    diversity_pr_local = "tables/diversity_pr_local.csv"
  params:
    diversity_path = "src/scripts/nucleotide_diversity.py",
    threads = 20
  shell: 
    """
    echo "Computing nucleotide diversity for cds_global_trimal..." 
    parallel -j {params.threads} 'python {params.diversity_path} -f {{}} -o {output.diversity_cds_global}' ::: {input.trimal_cds_global_align}/*.fna > cds_global.log 2> cds_global.err

    echo "Computing nucleotide diversity for promoter_global_trimal..."
    parallel -j {params.threads} 'python {params.diversity_path} -f {{}} -o {output.diversity_pr_global}' ::: {input.trimal_pr_global_align}/*.fna 

    echo "Computing nucleotide diversity for cds_local_trimal..." 
    parallel -j {params.threads} 'python {params.diversity_path} -f {{}} -o {output.diversity_cds_local}' ::: {input.trimal_cds_local_align}/*.fna 

    echo "Computing nucleotide diversity for promoter_local_trimal..." 
    find {input.trimal_pr_local_align} -name "*.fna" | parallel -j {params.threads} 'python {params.diversity_path} -f {{}} -o {output.diversity_pr_local}' 
    """

rule pangenome_occupancy_table:
  input:
    pangenome_matrix_shell= "out/clusters_cds/pangenome_matrix_t4__shell_list.txt",
    pangenome_matrix_softcore= "out/clusters_cds/pangenome_matrix_t4__softcore_list.txt"

  output:
    pangenome_table = "tables/pangenome.csv",
    occupancy_table = "tables/occupancy_number.csv"

  shell: 
    """
    #load pangenome table 
    sed 's/$/,shell/' {input.pangenome_matrix_shell} > tables/pangenome_matrix_t0__shell_list.txt
    sed 's/$/,core/' {input.pangenome_matrix_softcore} > tables/pangenome_matrix_t0__softcore_list.txt
    cat tables/pangenome_matrix_t0__shell_list.txt tables/pangenome_matrix_t0__softcore_list.txt > {output.pangenome_table}
    
    #load occupancy table 
    cat out/cds_est_homologues/BdistachyonABR2337v1_4taxa_algOMCL_e1_.cluster_list |\
    grep "taxa"|cut -d " " -f 4,6|perl -p -e 's/taxa=/,/'|\
    perl -lane 'print $F[1],$F[0]'> {output.occupancy_table}
    """

rule nucleotide_diversity_global_plot:
  input:
    script_nucleotide_diversity_plot = "src/scripts/nucleotide_diversity_plot.R"
  output:
    q_boxplot = "figures/nucleotide_diversity_global.png"
  
  shell:
    "Rscript {input.script_nucleotide_diversity_plot}  /home/ammar/ammar/snakemake_improve"


rule Mash_distance_estimation:
  input:
    cds_cluster= "out/clusters_cds/",
    promoter_cluster= "out/uncom_data/promoter_cluster/"
  output:
    cds_mash = directory("out/cds_mash/"),
    promoter_mash = directory("out/promoter_mash/"),
  params:
    script_mash = "src/scripts/mash2cluster.py",
    threads = 20

  shell:
    """
    mkdir -p {output.cds_mash} {output.promoter_mash}
    echo "Computing mash distance for cds..." 
    find {input.cds_cluster} -name "*.fna" | parallel -j {params.threads} 'python {params.script_mash} {{}} {output.cds_mash}/'

    echo "Computing mash distance for promotrs..." 
    find {input.promoter_cluster} -name "*.fna" | parallel -j {params.threads} 'python {params.script_mash} {{}} {output.promoter_mash}/'
    """ 


rule median_mash_distance:
  input:
    cds_mash = "out/cds_mash/",
    promoter_mash = "out/promoter_mash/",
  output:
    cds_median_mash = "tables/total_med_mash_dist_cds.csv",
    promoter_median_mash = "tables/total_med_mash_dist_pr.csv",
  params:
    script_median_mash = "src/scripts/median_mash_distance.py",
    threads = 20

  shell:
    """
    #how we extract mash distance from mash output
    for i in $(ls {input.cds_mash}/*_mash) ; do
      filename=$(basename "$i")
      cat $i | perl -pe 's/\t/,/g' | cut -d "," -f 3 | paste -sd , > {input.cds_mash}/$filename.dist
    done

    for i in $(ls {input.promoter_mash}/*_mash) ; do
      filename=$(basename "$i")
      cat $i | perl -pe 's/\t/,/g' | cut -d "," -f 3 | paste -sd , > {input.promoter_mash}/$filename.dist
    done

    #how we calculate the median script for all cluster and save it 
    find {input.cds_mash} -name "*.dist" | parallel -j {params.threads} 'python {params.script_median_mash} {{}} >> {output.cds_median_mash}'
    find {input.promoter_mash} -name "*.dist" | parallel -j {params.threads} 'python {params.script_median_mash} {{}} >> {output.promoter_median_mash}'
    """

rule mash_distance_plot:
  input:
    script_mash_plot = "src/scripts/mash_dist_plot.R"
  output:
    mash_boxplot = "figures/mash_distance.png",
    data_summary_mash = "tables/data_summary_mash.csv"
  
  shell:
    """
    echo "Running mash_distance_plot rule..."
    Rscript {input.script_mash_plot} /home/ammar/ammar/snakemake_improve
    rm  Rplots.pdf
    """

rule promoters_soft_extraction:
  input:
    script_soft_pr = "src/scripts/promoters_soft_extraction.sh"
  output:
    data_pr_soft = directory ("out/uncom_data/promoter_soft/")
  
  shell:
    "bash {input.script_soft_pr}"


rule promoters_soft_clustering:
  input:
    cds_ids = "out/uncom_data/cds/clusters_cds_ids/",
    script_seq_extraction = "src/scripts/seq_extraction_by_id.py",
    promoter_soft_dir = "out/uncom_data/promoter_soft/"
  output:
    promoter_soft_cluster = directory("out/uncom_data/promoter_soft_cluster/")
  shell:
    """
    mkdir -p {output.promoter_soft_cluster}/
    python {input.script_seq_extraction} -id {input.cds_ids} -f {input.promoter_soft_dir} -o {output.promoter_soft_cluster}/
    """



rule mask_cds_soft:
  input:
    cds_files = "out/uncom_data/cds/",
    script_mask_soft = "./src/scripts/masking_repeats/plant-scripts/repeats/Red2Ensembl.py"

  output:
    cds_masked = directory("out/uncom_data/cds_masked/"),

  shell:
    """
    # Run repeat masking and annotation tools for all fasta files
    for i in {input.cds_files}/*.fa; do
      NAME_base=$(basename "$i")
      echo "$NAME_base"
      {input.script_mask_soft} "$i" "out/uncom_data/cds_masked/${{NAME_base}}_file" \
      --msk_file "out/uncom_data/cds_masked/${{NAME_base}}.sm.fna" \
      --bed_file "out/uncom_data/cds_masked/${{NAME_base}}.bed" --cor 16 
    done
    """


rule cds_sm_clustering:
  input:
    cds_ids = "out/uncom_data/cds/clusters_cds_ids/",
    script_seq_extraction = "src/scripts/seq_extraction_by_id.py",
    cds_sm_dir = "out/uncom_data/cds_masked"
  output:
    cds_sm_cluster = directory("out/uncom_data/cds_cluster/")
  shell:
    """
    mkdir -p {output.cds_sm_cluster}/
    {input.script_seq_extraction} -id {input.cds_ids} -f {input.cds_sm_dir} -o {output.cds_sm_cluster}/
    """
  
rule count_TEs_soft_masked:
  input:
    promoter_soft_cluster = "out/uncom_data/promoter_soft_cluster/",
    cds_masked = "out/uncom_data/cds_cluster/",
    script_count_TEs = "src/scripts/count_masked_regions.py"

  output:
    count_masked_regions_promoter = "tables/count_te_sm_pr.csv",
    count_masked_regions_cds = "tables/count_te_sm_cds.csv"
  
  params:
    threads = 20

  shell:
    """
    find {input.promoter_soft_cluster} -name "*.fna" | \
    parallel -j {params.threads} 'python {input.script_count_TEs} {{}}\
    >> {output.count_masked_regions_promoter}'
    find {input.cds_masked} -name "*.fna" | \
    parallel -j {params.threads} 'python {input.script_count_TEs} {{}}\
    >> {output.count_masked_regions_cds}'
    """

rule TEs_content_plot:
  input:
    script_TEs_content_plot = "src/scripts/te_percentage_content_plot.R"
  output:
    TEs_content_boxplot = "figures/TEs_content.png"
  
  shell:
    """
    echo "Running TEs_content_plot rule..."
    Rscript {input.script_TEs_content_plot} /home/ammar/ammar/snakemake_improve
    rm  Rplots.pdf
    """

