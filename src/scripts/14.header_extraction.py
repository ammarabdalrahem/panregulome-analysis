#!/usr/bin/env python3
import argparse
#import glob


def run(args):
    inputfile = args.fasta_file
    outputfile = args.output_file
  
    fasta = open(inputfile, "r")
    output= open (outputfile,"w")
  

    headers = []

    for line in fasta:
      if line[0] == '>':
        header = line.split(" ")[0][1:]
        if '.' in header:
          header = header.split('.')[0]
        headers.append(header)
    for i in headers:
      output.write(i + "\n")

    output.close()

  
                    
def main():
    parser = argparse.ArgumentParser(description='print FASTA headers')
    parser.add_argument('-i', type=str, dest="fasta_file", required=True, help="insert name of fasta file")
    parser.add_argument('-o', type=str, dest="output_file", required=True, help="insert name of the output txt file")
    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)



if __name__ == "__main__":
    main()
