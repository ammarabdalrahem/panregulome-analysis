#! /usr/bin/bash    
    
    mkdir -p out/uncom_data/promoter_soft
    #have total softmasked genome sequences
    cp src/scripts/masking_repeats/plant-scripts/repeats/*.fa.sm.fna out/uncom_data/promoter_soft/

    #chnage name of softmasked genome sequences
    for file in out/uncom_data/promoter_soft/*.fa.sm.fna; do
      mv "$file" "${file%.fa.sm.fna}.fa"
    done

    #have total gff3 files
    cp out/uncom_data/genomes/*.gff3 out/uncom_data/promoter_soft/
  
  
    #extract (gene) from gff files |cut the gene ID |extract id,strand (+ or -),start and end coordinates
    # for + std sub 500 from sta but is sta is zero let it equal zero  ,and if - add 500 for end str

    ls -1 out/uncom_data/promoter_soft/*.gff3|
        sort -u |
        while read i
        do
            perl -lane 'if($F[2] eq "gene" && $F[8]=~/ID=([^;]+)/){ ($id,$sta,$end,$str)=($1,$F[3],$F[4],$F[6]);if($str eq "+"){$end=$sta; $sta-=500; 
            if($sta<0){ $sta=0 }} else { $sta=$end; $end+=500 } print "$F[0]\t$sta\t$end\t$id\t0\t$str" }' $i > $i.bed 
        done

    for i in $(ls out/uncom_data/promoter_soft/*gff3.bed | rev | cut -c9- |rev); do mv ${i}gff3.bed ${i}fa.bed;done


    #flank coordinates within (bedtools slop uning fasta file) to slove (the beyond the length) error make length not more than chr. size 
    ls -1 out/uncom_data/promoter_soft/*.fa|
        sort -u |
        while read i
        do

            perl -lne 'if (/^>(\S+)/){ $id=$1 } else { $len{$id} += length($_) } END{ foreach $id (sort keys(%len)){ print "$id\t$len{$id}" } }' $i > $i.chr.genome
            bedtools slop -l 0 -r 0 -s -i $i.bed -g $i.chr.genome > $i.flank.bed

        done

    #Run bedtools for promoter sequence extraction from bed file coordinates
    #input -fi fasta file -bed bed flank file, output -fo name new file

    ls -1 out/uncom_data/promoter_soft/*.fa|
        sort -u |
        while read i
        do

    bedtools getfasta -fi $i -bed $i.flank.bed -s -name -fo $i.pr.soft.fa

        done
