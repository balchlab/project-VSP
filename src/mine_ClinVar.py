
from __future__ import print_function
import time


def find_good_lines(text_file):

    i='0'
    for line in text_file:

        if line.startswith('#'):
            yield line
            continue

        elif 'indel' in line or 'deletion' in line:
            continue

        elif 'UniProtKB (protein)' in line:

            yield line



with open('../Data/ClinVAR_variant_summary.txt', 'rt') as tsvin, open('../Data/Clinvar_processed.txt', 'w') as txtout:
    Count_pathogenic = 0
    Count_Uncertain = 0
    Count_benign = 0
    for good_line in find_good_lines(tsvin):

        if 'Pathogenic' in good_line and 'Benign' not in good_line:
            Count_pathogenic +=1
        elif 'Uncertain significance' in good_line and 'Pathogenic' not in good_line:
            Count_Uncertain +=1
        elif 'Benign' in good_line and 'Pathogenic' not in good_line:
            Count_benign +=1


print ('Pathogenic:', Count_pathogenic)
print('Uncertain:', Count_Uncertain)
print('Benign:', Count_benign)

