
replacements = {
'Gly':'G',
'Ala':'A',
'Leu':'L',
'Met':'M',
'Phe':'F',
'Trp':'W',
'Lys':'K',
'Gln':'Q',
'Glu':'E',
'Ser':'S',
'Pro':'P',
'Val':'V',
'Ile':'I',
'Cys':'C',
'Tyr':'Y',
'His':'H',
'Arg':'R',
'Asn':'N',
'Asp':'D',
'Thr':'T'}



lines = []
with open('Mutant_list.txt') as infile:
    for line in infile:
        for src, target in replacements.items():
            line = line.replace(src, target)
        lines.append(line)
with open('Converted_Mutant_list.txt', 'w') as outfile:
    for line in lines:
        outfile.write(line)





# f1 = open('Mutant_list', 'r')
# f2 = open('Mutant_list.txt.tmp', 'w')
# for line in f1:
#     f2.write(line.replace('Ser', 'S'))
# f1.close()
# f2.close()