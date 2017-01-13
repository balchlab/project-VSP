
from __future__ import print_function
import time

import csv
import pandas as pd
import re


##Mine and clean up ClinVar txt file containing clinical information on pathogenic and neutral variants across genetic disease

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



def Convert_amino_acids(text):
    if type(text) == str:
        for i, j in replacements.items():
            text = text.replace(i, j)

        return (text)
    else:
        return ('nan')


def find_good_lines(text_file):

    for line in text_file:

        if line.startswith('#AlleleID'):
            print('found column names')
            yield line

        elif 'UniProtKB (protein)' in line:

            yield line

        elif 'indel' in line or 'deletion' in line:
            continue




def mine_clinvar():
    with open('../Data/ClinVAR_variant_summary.txt', 'rt') as tsvin, open('../Data/Clinvar_processed.csv', 'w') as csvout:
        write_csv = csv.writer(csvout)
        Count_pathogenic = 0
        Count_Uncertain = 0
        Count_benign = 0

        for good_line in find_good_lines(tsvin):


            if good_line.startswith('#AlleleID'):
                first_row = good_line.rstrip().split('\t')
                write_csv.writerows([first_row])


            elif 'Pathogenic' in good_line and 'Benign' not in good_line:

                Count_pathogenic +=1
                fields = good_line.rstrip().split('\t')
                write_csv.writerows([fields])


            elif 'Uncertain significance' in good_line and 'Pathogenic' not in good_line:

                Count_Uncertain +=1
                fields = good_line.rstrip().split('\t')
                write_csv.writerows([fields])

            elif 'Benign' in good_line and 'Pathogenic' not in good_line:

                Count_benign +=1
                fields = good_line.rstrip().split('\t')
                write_csv.writerows([fields])





    print ('Pathogenic:', Count_pathogenic)
    print('Uncertain:', Count_Uncertain)
    print('Benign:', Count_benign)

def build_clinvar_data_frame ():
    with open ('../Data/Clinvar_processed.csv', 'rt') as csvin, open('../Data/Clinvar_built.csv', 'w') as csvout:
        df = pd.read_csv(csvin)
        #print (df.head(5))


        sub_df = df[['Name', 'ClinicalSignificance', 'ClinSigSimple', 'GeneSymbol', 'OtherIDs']]
        #print (sub_df)


        sub_df['OtherIDs']= sub_df['OtherIDs'].str.replace('protein\):','#')

        sub_df['UniProt'] = sub_df['OtherIDs'].str.extract('#(.*)#')


        sub_df['Mutation'] = sub_df['Name'].str.replace('\(p.','\(#p.')

        sub_df['Mutation'] = sub_df['Mutation'].str.extract('#(.*)\)')

        sub_df['Mutation_2'] = sub_df['Mutation'].apply(Convert_amino_acids)

        sub_df = sub_df[sub_df.Mutation_2 !='nan']

        sub_df['Mutation_3'] = sub_df['Mutation_2'].str.replace('p.','')

        sub_df['AA_position'] = sub_df['Mutation'].str.extract('(\d+)').astype(int)

        sub_df = sub_df.drop_duplicates('Mutation', take_last=True)

        sub_df.reset_index(drop=True, inplace = True)

        print (sub_df.head(5))
        sub_df.to_csv(csvout)

def main():
    mine_clinvar()
    build_clinvar_data_frame()


if __name__ == '__main__':
    main()
