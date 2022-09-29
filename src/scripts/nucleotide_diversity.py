#!/usr/bin/env python3


import glob
import argparse

def run(args):
    fasta_file=args.fasta_file
    outputfile = args.output_directory


    path=fasta_file
    output_path=outputfile

    # extraction sequence in fasta format and save it in dictionary
    seqs={}
    for filename in glob.glob(path, recursive=True):
        with open(filename) as file:
            file_name=str(filename).split("/")[-1]
            for line in file :
                if line[0] == ">":                      #read the header
                    header = line.strip()                      #header = gene id
                    seqs[header] = ""                   #take gene id as main key
                else:
                    seqs[header] += line[:-1]           #add each line after header as value fo header key




    sequance = seqs.values()                            #  Retrieve dictionary values
    MSA = list(sequance)
    MSAlength=len(MSA[1])                               # assumes all sequences in MSA have same length
    sequences_number=len(MSA)                                # number of sequences(taxa) in the cluster
    all_pi = []
    bias_estimator= (sequences_number/(sequences_number-1))

    for seqi in range(len(MSA)):                        # pairwise comparisons
        for seqj in range(seqi+1,len(MSA)):
            different=0
            totalaln=0
            for bp in range(MSAlength):
                if MSA[seqi][bp] in ['n', '-'] or MSA[seqj][bp] in ['n', '-']: continue #avoid gaps and Ns
                totalaln = totalaln+1

                if MSA[seqi][bp] != MSA[seqj][bp]:
                    different=different+1 #number of different sites

            #in case sequence contain all gaps or Ns to avoid error assumption totalaln equal 1
            if totalaln==0 :
                 totalaln=1
            #frequency of each sequence in population qual 1/number of total sequance
            frequency=(1/sequences_number)
            #pi is the nucleotide diversity for any pair of sequences
            pi=((different/totalaln)*frequency*frequency)
            all_pi.append (pi)
    #from Nucleotide diveristy Equation
    nucletodie_diversity=bias_estimator*(sum(all_pi))
    print(file_name,",",nucletodie_diversity, file=open(output_path, "a"))



def main():
    parser = argparse.ArgumentParser(description='calculate nucleotide diversity')
    parser.add_argument('-f', type=str, dest="fasta_file", required=True, help="insert name of the output txt file")
    parser.add_argument('-o', type=str, dest="output_directory", required=True, help="insert name of the output txt file")
    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)



if __name__ == "__main__":
    main()
