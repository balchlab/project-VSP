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
GENE = "PKD2"
FILENAME = "PKD2PROVEANScores.csv"

# change directory to working with DAta
os.chdir("../Data/")
cwd = os.getcwd()
#df = pd.read_csv('foo.csv', index_col=0)

def findPROVEANscores(protein_ID):
    #read from tsv.gz file
    with gzip.open('PROVEAN_scores_ensembl66_human.tsv.gz','rt') as tsvin, open(FILENAME, 'wt') as csvout:
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




def combineExACtoPROVEAN(input):


    df = (pd.read_csv(input))

    short_df = df.drop(df.columns[[0,1]], axis=1)
    short_df['AA'] =(short_df.T.idxmax())
    short_df['position'] =df['position'].apply(str)
    short_df['Mutation'] = short_df['AA'].astype(str) + short_df['position'].astype(str)
    df.to_csv(FILENAME, sep='\t')

    print(short_df)





def main ():

    findPROVEANscores(ENSP)

    combineExACtoPROVEAN(FILENAME)

if __name__ == '__main__':
    main()




