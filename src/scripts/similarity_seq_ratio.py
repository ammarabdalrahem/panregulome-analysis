#!/usr/bin/env python3

#excluding poorly aligned sequence after global alignment

import glob
import argparse
import statistics

def run(args):
    fasta_file=args.fasta_file  # define fasta file name
    outputfile = args.output_directory # define the name output fasta file

    path=fasta_file
    output_path=outputfile


    print("reading sequences...")
    seqs={} #creat dictionary to all fasta files with seq. header as key and sequance as value.
    headers_list= []
    for filename in glob.glob(path, recursive=True):
        file_name=str(filename).split("/")[-1] # file name
        with open(filename) as file:
            for line in file :
                if line[0] == ">": #read the header
                    header = line.strip() #header = gene id
                    headers_list.append (header) # creat list of sequance's headers with the same order in fasta file.
                    seqs[header] = "" #take gene id as main key
                else:
                    seqs[header] += line[:-1] #add each line after header as value fo header key


    MSAlength=len (seqs[headers_list[0]]) #assumes all sequences in sequance have same length
    ratio_dic={}

    for seqi in range(len(headers_list)):# calculate the identity percent
        for seqj in range(seqi+1,len(headers_list)): # to compare all possibilities
            identities=0
            totalaln=0
            for bp in range(MSAlength):
                if seqs[ headers_list[seqi] ][bp] in ['n', '-'] or seqs [headers_list [seqj] ][bp] in ['n', '-']: continue
                totalaln = totalaln+1
                if seqs [headers_list [seqi]] [bp] == seqs [headers_list[seqj]][bp]:
                    identities=identities+1

            if totalaln == 0 : # in case no match, to avoid the error.
                totalaln=1


            ratio=((identities/totalaln)*100) #calculate the ratio
            ratio_dic[ headers_list[seqi],headers_list [seqj]] = ratio #save ratio in dictinary i,j & j,i
            ratio_dic[ headers_list [seqj],headers_list [seqi]] = ratio


    f_seqs = {k: list(v.split(" ")) for k,v in seqs.items()} # create final dic to remove poor alignment sequances.


    # calculate the median of identities
    for seqi in range(len(headers_list)):
        ratios = []
        for seqj in range(len(headers_list)):
            if( headers_list [seqi] ==  headers_list [seqj]): continue #avoid duplication
            ratios.append ( ratio_dic[ headers_list[seqi],headers_list [seqj]  ])
            ratios.append ( ratio_dic[ headers_list[seqj],headers_list [seqi]  ])

            med = statistics.median(ratios) # count the median
        print(med,file=open(output_path, "a"))
    print("done")




def main():
    parser = argparse.ArgumentParser(description='filtere the poor alignment sequances')
    parser.add_argument('-f', type=str, dest="fasta_file", required=True, help="insert name of the output fna file")
    parser.add_argument('-o', type=str, dest="output_directory", required=True, help="insert name of the output fna file")
    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)



if __name__ == "__main__":
    main()
