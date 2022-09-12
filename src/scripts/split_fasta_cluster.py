#!/usr/bin/env python3

#excluding poorly aligned sequence after global alignment

import glob
import sys
import os
import shutil



input_file=sys.argv[1]
output_dr = str (input_file) + "_temp/"
os.mkdir (output_dr) 
print(output_dr, "created")


print("reading sequences...")
seqs={} #creat dictionary to all fasta files with seq. header as key and sequance as value.
headers_list= []
tmpids=[]
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
	tempfilename= output_dr + tmpfilename
	tmpids.append(tempfilename)


for i in range(len(headers_list)):# calculate the identity percent	
	for k in seqs:
		print(k,"\n",seqs[k],file=open(tmpids[i] , "a"))

"""
for seqi in range(len(headers_list)):# calculate the identity percent
	for seqj in range(seqi+1,len(headers_list)): # to compare all possibilities
		print (seqi,"temporary",input_file,sep='')
		print (seqj,"temporary",input_file,sep='')

	
	

if  output_dr.endswith('_temp'):
	shutil.rmtree(output_dr)
	print ("temporary directory deleted")
"""
