#! /usr/bin/bash

# install bedtools (v2.29.1)
cd ./../scripts
wget https://github.com/arq5x/bedtools2/releases/download/v2.29.1/bedtools-2.29.1.tar.gz
tar -zxvf bedtools-2.29.1.tar.gz
cd bedtools2
make
cd ./../../work

#convert gff file to bed file 

#extract (gene) from gff files |cut the gene ID |extract id,strand (+ or -),start and end coordinates
# for + std sub 500 from sta but is sta is zero let it equal zero  ,and if - add 500 for end str

ls -1 *.gene.gff3|
    sort -u |
    while read i
    do
        perl -lane 'if($F[2] eq "gene" && $F[8]=~/ID=([^;]+)/){ ($id,$sta,$end,$str)=($1,$F[3],$F[4],$F[6]);
            if($str eq "+"){ $sta-=500; if($sta<0){ $sta=0 } } else { $end+=500 }  
            print "$F[0]\t$sta\t$end\t$id\t0\t$str" } ' $i> $i.bed 
    done 


#change the extention name to be as fasta file name plus .bed for easy run in bedtools

for f in ./*.bed; do      mv -- "$f" "${f%.*.1.gene.gff3.bed}.fa.bed"; done


#flank coordinates within (bedtools slop uning fasta file) to slove *the beyond the length)error make length not more than chr. size 
ls -1 *.fa|
    sort -u |
    while read i
    do

	perl -lne 'if(/^>(\S+)/){ $id=$1 } else { $len{$id} += length($_) } END{ foreach $id (sort keys(%len)){ print "$id\t$len{$id}" } }' $i> $i.chr.genome
	bedtools slop -l 0 -r 0 -s -i $i.bed -g $i.chr.genome > $i.flank.bed
done


#Run bedtools for promoter sequence extraction from bed file coordinates
#input -fi fasta file -bed bed flank file, output -fo name new file
ls -1 *.fa|
    sort -u |
    while read i
    do

bedtools getfasta -fi $i -bed $i.flank.bed -s -fo $i.promotor.out

done
