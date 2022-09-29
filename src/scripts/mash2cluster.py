#!/usr/bin/env python3

# TODO

import glob
import sys
import os
import shutil
import errno
import subprocess

mash_exe = '/home/contrera/soft/mash-Linux64-v2.3/mash'
if not (os.path.isfile(mash_exe)):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), mash_exe)

# TODO check whether mash is available

input_file=sys.argv[1] #input file name

# TODO: exit with example call if no input is passed

# work our output directory
output_dr = str (input_file) + "_temp/" #the name of directory with add temp to delete it at the end

os.mkdir (output_dr)  #creat directory as name of original file + temp
print(output_dr, "created")


print("reading sequences...")
seqs={} #creat dictionary to all fasta files with seq. header as key and sequance as value.
headers_list= [] #list of headers to keep it arranget
tmpids=[] #list with new files names
for filename in glob.glob(input_file, recursive=True):
        file_name=str(filename).split("/")[-1] # file name
        with open(filename) as file:
                for line in file :
                        if line[0] == ">": #read the header
                                header = line.strip() #header = gene id
                                headers_list.append (header) # creat list of sequance's headers with the same order in fasta file.
                                seqs[header] = "" #take gene id as main key
                        else:
                                seqs[header] += line[:-1] #add each line after header as value fo header key


for i in range (len(headers_list)):
        tmpfilename = headers_list[i]
        tmpfilename=str(headers_list[i]).split(">")[-1] # file name
        tmpfilename= tmpfilename.split(" ")[0]
        tmpfilename= tmpfilename+ ".fna"
        #tempfilename= output_dr + tmpfilename
        tmpids.append(tmpfilename)


for i in range(len(headers_list)):# calculate the identity percent
        ID=headers_list[i]
        print(ID,"\n",seqs[ID],file=open(output_dr + tmpids[i] , "a"))


for seqi in range(len(headers_list)):# calculate the Kmer distance
        for seqj in range(seqi+1,len(headers_list)): # to compare all possibilities

                file1=tmpids[seqi]
                file2=tmpids[seqj]
                print (file1,file2)

                # open new log file
                logfile = sys.argv[2] + str (input_file) + "_mash"
                try:
                        logfile = open(logfile,"a")
                except OSError as error:
                        print("# ERROR: cannot create file ", log_filepath, error)

                # put together call to mash
                cmd = (
                        mash_exe
                        + " dist "
                        +" "
                        + output_dr
                        + str(file1)
                        + " "
                        + output_dr
                        + str(file2)
                )

                # run mash and capture stdout


                try:
                        print("# mash command: ", cmd)
                        osresponse = subprocess.check_call(cmd.split(), stdout=logfile)
                except subprocess.CalledProcessError as err:
                        print("# ERROR: cannot run mash ", err.returncode)
                        logfile.close()
                        #return rpt_files
                finally:
                        logfile.close()


if  output_dr.endswith('_temp/'):
        shutil.rmtree(output_dr)
        print ("temporary directory deleted")
