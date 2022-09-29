#!/usr/bin/env python3
#@AbdalrahemAmmar


#count soft-masked regions in sequance

import os
import sys
from Bio import SeqIO



input_file = open(sys.argv[1], 'r')
count=[]
total_length_seq=[]

for record in SeqIO.parse(input_file, "fasta") :
    length = len(record.seq)
    total_length_seq.append (length)
    a_count = record.seq.count('a')
    c_count = record.seq.count('c')
    g_count = record.seq.count('g')
    t_count = record.seq.count('t')
    total=a_count+c_count+g_count+t_count
    count.append (total)

total_seq = sum (total_length_seq) #total sequance length
masked_count = sum(count) #total lowercase (soft masked)
precentage_te = int((masked_count/total_seq)*100) #calculate the % of soft-masked regions


print (total_seq,",",masked_count,",",precentage_te)
