from __future__ import print_function
from chardet.universaldetector import UniversalDetector
import sys, traceback,subprocess,gzip,glob, tarfile,os, signal
import pandas as pd
import csv
import numpy as np
import numpy as np
from pandas import HDFStore,DataFrame
import h5py


# create (or open) an hdf5 file and opens in append mode


os.chdir("../Data/")
cwd = os.getcwd()

# def test_encoding(file_name):
#     detector = UniversalDetector()
#     with open(file_name, 'rb') as f:
#         for line in f:
#             detector.feed(line)
#             if detector.done:
#                  break
#         detector.close()
#     r = detector.result
#     return "Detected encoding %s with confidence %s" % (r['encoding'], r['confidence'])

with gzip.open('PROVEAN_scores_ensembl66_human.tsv.gz','rt') as tsvin, open('new.csv', 'wt') as csvout:
    tsvin = csv.reader(tsvin, delimiter='\t',quoting=csv.QUOTE_NONE)
    csvout = csv.writer(csvout)

    for row in tsvin:
        count = row[0]
        if count ==  "ENSP00000237596":
            csvout.writerows([row[0:8]])

#datatable=csvreader('PROVEAN_scores_ensembl66_human.tsv.gz')

#print (datatable)


#PRODATA = pd.read_csv('PROVEAN_scores_ensembl66_human.tsv.gz',delimiter="\t", quoting=csv.QUOTE_NONE, encoding='utf-8-sig', nrows = 100000)
#hdf.put('d1', PRODATA, format='table', data_columns=True)
#PRODATA.to_hdf('PROVEAN_scores.hdf5', '/data', append=True)


#store = HDFStore('PROVEAN.h5')
#store['df'] = PRODATA  # save it
#store['df']  # load it
#with h5py.File('PROVEAN.h5', 'w') as hf:

    #hf.create_dataset('d1', data = PRODATA, compression="lzf", compression_opts=9)


#print (store.df.head())