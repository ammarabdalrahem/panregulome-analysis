
import argparse

def run(args):
    ids_file = args.ids_file
    fasta_file=args.fasta_file
    outputfile = args.output_directory

#all input only dic paths
    path=fasta_file
    path_2=ids_file
    output_path=outputfile

    print("reading promoter sequences...")
    seqs={}
    for filename in glob.glob(path + '/*', recursive=True):
        with open(filename) as file:
            for line in file :
                if line[0] == ">":                      #read the header
                    header1 = line.split(">") [1]      #header = gene id
                    header = header1.strip()            # to remve new line for cds headers
                    sp_info=line.strip()                #species name
                    seqs[header] = {}                   #take gene id as main key
                    seqs[header][sp_info]= ""           # add seq. to second key (sp_info)

                else:
                    seqs[header][sp_info] += line[:-1]    #add each line after header as value fo header key

    print("done");


    # iterate through gene clusters, each a FASTA file
    for filename in glob.glob(path_2+ '/*', recursive=True):
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


def main():
    parser = argparse.ArgumentParser(description='print FASTA headers')
    parser.add_argument('-id', type=str, dest="ids_file", required=True, help="insert name of fasta file")
    parser.add_argument('-f', type=str, dest="fasta_file", required=True, help="insert name of the output txt file")
    parser.add_argument('-o', type=str, dest="output_directory", required=True, help="insert name of the output txt file")
    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)



if __name__ == "__main__":
    main()

# to run it for promoter    
# ./src/scripts/seq_extraction_by_id.py -f out/cds/cds_sm  -id out/cds/clusters_cds_ids -o out/cds/clusters_cds_sm/
# to change the name
# for f in ./*.txt; do  mv -- "$f" "${f%.fna.txt}.fna"; done
