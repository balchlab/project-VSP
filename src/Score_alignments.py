import sys

import collections, time
import pandas as pd


 # open and read fasta files
 # sperate names and store separately
 #

df = pd.read_csv("../Alignments/Km_data.csv")

Rubisco_data = df[['RbcL Sequence', 'Km CO2', 'Km O2']]
print (Rubisco_data.head(5))
Rubisco = {}
Rubisco = Rubisco_data.set_index('RbcL Sequence').T.to_dict('list')



fasta = {}
with open("../Alignments/Rubisco.txt") as file_one:
    for line in file_one:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            active_sequence_name = line[1:]
            if active_sequence_name not in fasta:
                fasta[active_sequence_name] = []
            continue
        sequence = line
        fasta[active_sequence_name].append(list(sequence))

i = 0

while i < 100: # this needs to be changed to total length of alignment

    str_list = []
    s = ""
    low_freq_AAs = {}
    # first calculate frequencies of each amino acid in position
    # then iterate through these positions again
    # if found match and freq less than 50%
    # print out stats and name of key
    
    for key in fasta:

        for x in (fasta[key]):

            s+=x[i-1]

    #print (s)
    Count = collections.Counter(s)

    for key, value in Count.items ():

        #print (i+1, key, value/35*100,'%')

        if 50 > value/35*100:
            #print ("low")
            #print ("key:", key)

            low_freq_AAs[key] = value/35*100

    #print (low_freq_AAs)


    for key in fasta:

        for x in (fasta[key]):

            if x[i-1] in low_freq_AAs:
                name = key.partition('|')[-1].rpartition('|')[0]
                Kms_list = []
                Kms_list = Rubisco[name]
                Kms = ' '.join(str(v) for v in Kms_list)
                print (i+1, x[i-1], low_freq_AAs[x[i-1]], name, Kms)


    i+=1

