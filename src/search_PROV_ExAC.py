from __future__ import print_function
from chardet.universaldetector import UniversalDetector
import sys, traceback,subprocess,gzip,glob, tarfile,os, signal
import pandas as pd
import csv
import numpy as np
import numpy as np
from pandas import HDFStore,DataFrame
import h5py



#Sets protein ID to search in dataframe
ENSP = "ENSP00000237596"
ENSG = "ENSG00000118762"
ENST = "ENST00000237596"
GENE = "PKD2"
FILENAME1 = "PKD2PROVEANScores.csv"
FILENAME2 = "PKD2ExACScores.csv"
FILENAME3 = "PKD2MutPredScores.csv"
UniProt = "Q13563"

# change directory to working with DAta
os.chdir("../Data/")
cwd = os.getcwd()
#df = pd.read_csv('foo.csv', index_col=0)

def findPROVEANscores(protein_ID):
    #read from tsv.gz file
    with gzip.open('PROVEAN_scores_ensembl66_human.tsv.gz','rt') as tsvin, open(FILENAME1, 'wt') as csvout:
        csvout = csv.writer(csvout)
        tsvin = csv.reader(tsvin, delimiter='\t',quoting=csv.QUOTE_NONE)
        for i in range(1):

            row1 = next(tsvin)
            print(row1)
            csvout.writerows([row1])
            i=+1



        for row in tsvin:

            count = row[0]
            if count == protein_ID:
                csvout.writerows([row[0:23]])




def formatPROVEAN(input):


    df = (pd.read_csv(input))

    short_df = df.drop(df.columns[[0,1]], axis=1)
    short_df['AA'] =(short_df.T.idxmax())
    short_df['position'] =df['position'].apply(str)
    short_df['Mutation'] = short_df['AA'].astype(str) + short_df['position'].astype(str)
    df.to_csv(FILENAME, sep='\t')

    print(short_df)

def mineExAC(SYMBOL):

    #read from tsv.gz file
    with gzip.open('ExAC.r0.3.1.sites.vep.table.gz','rt') as tsvin, open(FILENAME2, 'wt') as csvout:
        csvout = csv.writer(csvout)
        tsvin = csv.reader(tsvin, delimiter='\t',quoting=csv.QUOTE_NONE)
        print ('looking for this query: ',SYMBOL)
        for i in range(1):
            row1 = next(tsvin)
            print(row1)
            print('found ', len(row1), 'rows' )

            csvout.writerows([row1])
            i=+1


        variants = 0
        for row in tsvin:

            count = row[62] #row 60 is ENSG, 62 is SYMBOL, 61 is Feature


            if count == SYMBOL:
                variants +=1
                print(variants,' writing', SYMBOL, 'variant found in chromosome ', row[0])
                csvout.writerows([row[0:len(row1)]])

        print ('found ',variants, 'variants')

def mineMutPred(S):

    #read from tsv.gz file
    with gzip.open('MutPred.txt.gz','rt') as tsvin, open(FILENAME3, 'wt') as csvout:

        csvout = csv.writer(csvout)
        tsvin = csv.reader(tsvin, delimiter='\t',quoting=csv.QUOTE_NONE)
        print ('looking for this query: ',S)
        for i in range(1):
            row1 = next(tsvin)
            print(row1)
            print('found ', len(row1), 'rows' )

            csvout.writerows([row1])
            i=+1


        variants = 0
        for row in tsvin:


            count = row[1]
            #print (row[1],S)
            if count == S:
                variants +=1
                print(variants,' writing', S, 'variant scores ', row[0])
                csvout.writerows([row[0:len(row1)]])

        print ('found ',variants , 'variant MutPred Scores')


def main ():

    #findPROVEANscores(ENSP)

    #formatPROVEAN(FILENAME)
    #mineExAC(GENE)
    mineMutPred(UniProt)

if __name__ == '__main__':
    main()




