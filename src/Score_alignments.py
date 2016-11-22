import sys
import collections



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
while i < 100:

    str_list = []
    s = ""
    for key in fasta:
        for x in (fasta[key]):
            s+=x[i]
    Count = collections.Counter(s)
    for key, value in Count.items ():
        print (i+1, key, value/35*100,'%')

    i+=1




