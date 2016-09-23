from __future__ import print_function
from chardet.universaldetector import UniversalDetector
import sys, traceback,subprocess,gzip,glob, tarfile,os, signal
import pandas as pd
import csv
import numpy as np
import numpy as np
from pandas import HDFStore,DataFrame
import h5py
import zipfile
import odo

#Sets protein ID to search in dataframe
ENSP = "ENSP00000237596"
ENSG = "ENSG00000186868"
ENST = "ENST00000237596"
GENE = "PKD2"
FILENAME1 = "PKD2PROVEANScores.csv"
FILENAME2 = "PKD2ExACScores.csv"
FILENAME3 = "PKD2MutPredScores.csv"
FILENAME4 = "dbNSFP_output.csv"
UniProt = "Q13563"
Chr = "chrM"

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
def mine_dbNSFP(Chr, ENSG):
    ChrFilesDict = {
        'chr1':'dbNSFP3.2c_variant.chr1',
        'chr2':'dbNSFP3.2c_variant.chr2',
        'chr3':'dbNSFP3.2c_variant.chr3',
        'chr4':'dbNSFP3.2c_variant.chr4', #ENSG is row 19"
        'chr5':'dbNSFP3.2c_variant.chr5',
        'chr6':'dbNSFP3.2c_variant.chr6',
        'chr7':'dbNSFP3.2c_variant.chr7',
        'chr8':'dbNSFP3.2c_variant.chr8',
        'chr9':'dbNSFP3.2c_variant.chr9',
        'chr10':'dbNSFP3.2c_variant.chr10',
        'chr11':'dbNSFP3.2c_variant.chr11',
        'chr12':'dbNSFP3.2c_variant.chr12',
        'chr13':'dbNSFP3.2c_variant.chr13',
        'chr14':'dbNSFP3.2c_variant.chr14',
        'chr15':'dbNSFP3.2c_variant.chr15',
        'chr16':'dbNSFP3.2c_variant.chr16',
        'chr17':'dbNSFP3.2c_variant.chr17',
        'chr18':'dbNSFP3.2c_variant.chr18',
        'chr19':'dbNSFP3.2c_variant.chr19',
        'chr20':'dbNSFP3.2c_variant.chr20',
        'chr21':'dbNSFP3.2c_variant.chr21',
        'chrX':'dbNSFP3.2c_variant.chrX',
        'chrM':'dbNSFP3.2c_variant.chrM',
        }
    #read from tsv.gz file
    with zipfile.ZipFile('dbNSFPv3.2c.zip','r') as tsvin, open(FILENAME4, 'wt') as csvout:

        print (ChrFilesDict[Chr])
        #FileNames = tsvin.namelist()
        #print (len(FileNames))
        #print(FileNames)
        tp = pd.read_csv(tsvin.open(ChrFilesDict[Chr]), delimiter='\t' , quoting=csv.QUOTE_NONE)# iterator = True,  dtype=object)

        def iterrows (tp, chunksize):
            n = len(tp.index)
            i = 0
            while (i<n):
                slice = tp[i: i+chunksize]
                for row in slice:
                    yield(row)
                i+= chunksize
        colOne=tp.columns.tolist()
        print (colOne)
        for row in iterrows (tp, 10):
             print (row[0:10])
             count = row[0]
             print (row[0],ENSG)
             if count == ENSG:
                variants +=1
                print(variants,' writing', ENSG, 'variant scores ', row[0])
                csvout.writerows([row[0:len(row1)]])
        print (variants)






        #
        # for chunk in odo(tp, chunk(pd.DataFrame), chunksize = 100): print (chunk)
        #variants = 0
        #for chunk in tp:
            #print (chunk[0])
            #Data =

        #     for row in chunk:
        #         print (row[19])
        #         count = row[19]
        #         print (row[1],S)
        #         if count == ENSG:
        #             variants +=1
        #             print(variants,' writing', ENSG, 'variant scores ', row[0])
        #             csvout.writerows([row[0:len(row1)]])
        #
        # print (variants)
        # csvout = csv.writer(csvout)



def main ():

    #findPROVEANscores(ENSP)

    #formatPROVEAN(FILENAME)
    #mineExAC(GENE)
    #mineMutPred(UniProt)
    mine_dbNSFP(Chr, ENSG)

if __name__ == '__main__':
    main()




