#!/usr/bin/env python3

import glob
output_path="/home/ammar/ammar/brachy_snakemake/out/promoter_cluster/"

print("reading promoter sequences...")
seqs={}
for filename in glob.glob("/home/ammar/ammar/brachy_snakemake/out/proximal_promoter/*.fa", recursive=True):
    with open(filename) as file:
        for line in file :
            if line[0] == ">":                      #read the header
                header = line.split(".")[0][1:]     #header = gene id 
                sp_info=line.strip().split(":")[-4] #species name 
                seqs[header] = {}                   #take gene id as main key
                seqs[header][sp_info]= ""           # add seq. to second key (sp_info)

            else:
                seqs[header][sp_info] += line[:-1]    #add each line after header as value fo header key
print("done");

# iterate through gene clusters, each a FASTA file
for filename in glob.glob("/home/ammar/ammar/brachy_snakemake/out/cluster_ids/*.txt", recursive=True):
    with open(filename) as header_file:
        file_name=str(header_file).split("\'")[1].split("/")[-1]
        print("parsing " + file_name)

        # extract gene ids and species for this cluster
        headers = []
        for i in header_file:
            # headers in FASTA are gene ids such as Brdisv1ABR21003818m
            headers.append(i[:-1]) 
        #promoter cluster as gene cluster by gene ids
        output = open(output_path+file_name, 'w')
        for id in headers:
            if id in seqs:
                for sp, seq in seqs[id].items():
                    output.write(str(sp) +'\n'+ str(seq)+'\n')            

        output.close()
