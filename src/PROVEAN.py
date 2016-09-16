'''
Converts ExAC missense protein notation to single letter notation
TODO: get PROVEAN.

input: .txt file with "Arg13Pro" style notation
output: .txt file with  "R13P" notation
'''
#Dictionary with amino acids
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



def AAconverter(Input):
#returns converted text file list

    lines = []
    for line in Input:
        for src, target in replacements.items():
            line = line.replace(src, target)
        lines.append(line)

    with open('Converted_Mutant_list.txt', 'w') as outfile:
        for line in lines:
            outfile.write(line)
def main ():
    with open('Mutant_list.txt') as Input:
        AAconverter(Input)



if __name__ == '__main__':
    main()






# f1 = open('Mutant_list', 'r')
# f2 = open('Mutant_list.txt.tmp', 'w')
# for line in f1:
#     f2.write(line.replace('Ser', 'S'))
# f1.close()
# f2.close()