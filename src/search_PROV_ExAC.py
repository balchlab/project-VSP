'''
Missense variant miner:
find missense variants in ExAC, score them using PROVEAN, MutPred and dbNSFP
return csv file.



'''
from __future__ import print_function
from chardet.universaldetector import UniversalDetector
import sys, traceback,subprocess,gzip,glob, tarfile,os, signal
import pandas as pd
import csv
import numpy as np
from pandas import HDFStore,DataFrame
import h5py
import zipfile
import odo
from collections import OrderedDict
import collections
from time import sleep
import string
import json



VCF_HEADER = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']


#Sets protein ID to search in dataframe
ENSP = "ENSP00000003084"
ENSG = "ENSG00000001626"
#ENSG = "ENSG00000186868" #MAPT
#ENSG = "ENSG00000272636" #Diagnostic - beginning of Chr17
ENST = "ENST00000003084"

GENE = "CFTR"
FILENAME = "CFTR_PROV_extract.csv"
FILENAME1 = "CFTR_PROVEANScores.csv"
FILENAME2 = "CFTR_ExACScores.csv"
FILENAME3 = "CFTR_MutPredScores.csv"
FILENAME4 = "dbNSFP_output.csv"
FILENAME5 = "dbNSFP_extract.csv"
UniProt = "P13569"
Chr = "7"

# change directory to working with DAta
os.chdir("../Data/")
cwd = os.getcwd()
#df = pd.read_csv('foo.csv', index_col=0)
def findPROVEANscores(ENSP):
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
            #print (count, ENSP)
            if count == ENSP:
                csvout.writerows([row[0:23]])

def formatPROVEAN(input):
    df = (pd.read_csv(input))
    short_df = df.drop(df.columns[[0,1]], axis=1)
    short_df['AA'] =(short_df.T.idxmax())
    short_df['position'] =df['position'].apply(str)
    short_df['Mutation'] = short_df['AA'].astype(str) + short_df['position'].astype(str)
    df.to_csv(FILENAME, sep='\t')
    print(short_df)

def dataframeExAC(SYMBOL):
    print("Open")

    #read from tsv.gz file
    with gzip.open('ExAC.r0.3.1.sites.vep.vcf.gz','rt') as tsvin, open(FILENAME2, 'wt') as csvout:
        print("open")
        comments = count_comments(tsvin)
        result = OrderedDict()
        csvout = csv.writer(csvout)
        tsvin = pd.read_table(tsvin,  skiprows=comments,  names=VCF_HEADER, usecols=range(8))  #delimiter='|',quoting=csv.QUOTE_NONE, skiprows = comments)
        print ('looking for this query: ',SYMBOL)

        for i, line in enumerate(lines(tsvin)):
            for key in line.keys():
                if key not in result:
                    result[key] = [None] * i

            for key in result.keys():
                result[key].append(line.get(key, None))

        return pd.DataFrame(result)

def generatorExAC (filename, Chr):
    print('searching Chrom:', Chr, 'in file: ',filename)
    lines(filename, Chr)




        # for i in range (1):
        #     row1 = next(tsvin)
        #     print (row1)
        #     print ('found', len(row1), 'rows')
        #     csvout.writerows([row1])
        #     i+=1



        # for line in tsvin:
        #
        #     if line.startswith('#'):
        #         continue
        #     else:
        #         parse(line)

def mineMutPred(ENST, UniProt):
    #TODO: Have to search both ENST and UniPROT codes in same query
    #read from tsv.gz file
    with gzip.open('MutPred.txt.gz','rt') as tsvin, open(FILENAME3, 'wt') as csvout:

        csvout = csv.writer(csvout)
        tsvin = csv.reader(tsvin, delimiter='\t',quoting=csv.QUOTE_NONE)
        print ('looking for this query: ',ENST, UniProt)
        for i in range(1):
            row1 = next(tsvin)
            print(row1)
            print('found ', len(row1), 'rows' )
            csvout.writerows([row1])
            i=+1
        variants = 0
        for row in tsvin:
            count = row[1]
            #print(count,ENST, UniProt)
            if count == ENST or count == UniProt:
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
        tp = pd.read_csv(tsvin.open(ChrFilesDict[Chr]), delimiter='\t' , quoting=csv.QUOTE_NONE, iterator = True,  dtype=object, chunksize = 40)# header = None)
        writer = csv.writer(csvout)
        for i in range(1):
            row1 = next(tp)
            print(row1)
            print('found ', len(row1), 'rows' )
            writer.writerows([row1[1:len(row1)]])
            i=+1
        print (ChrFilesDict[Chr])

        variants = 0
        for chunk in tp:
            row = next(chunk.itertuples())
            count = row[20]
            #print (row[1:len(row)])
            if count == ENSG:
                variants +=1
                print(variants,' writing', ENSG, 'variant scores ', row[0])
                writer.writerows([row[1:len(row)]])
        print (variants)

def extract_dbNSFP(file):

    with open(file,'rt') as tsvin, open(FILENAME5, 'wt') as csvout:
        dict = []
        df = pd.read_csv(tsvin, delimiter ='\t')
        print (df)
        for i in range(5):
            row = next(df)
            print (row)
            dict.append(row)
            i+=1
        print(dict)

def count_comments(filename):
    """Count comment lines (those that start with "#") in an optionally
    gzipped file.
    :param filename:  An optionally gzipped file.
    https://gist.github.com/slowkow/6215557
    """
    comments = 0
    print("comm")
    #fn_open = gzip.open if filename.endswith('.gz') else open
    #with fn_open(filename) as fh:
    for line in filename:
        if line.startswith('#'):
            comments += 1
            print(line)
        else:
            print ("no more")
            break
    print (comments)
    return comments

def parse(line):
    """Parse a single VCF line and return an OrderedDict.
    https://gist.github.com/slowkow/6215557
    """
    result = OrderedDict()

    fields = line.rstrip().split('\t')

    # Read the values in the first seven columns.
    for i, col in enumerate(VCF_HEADER[:7]):
        result[col] = _get_value(fields[i])

    # INFO field consists of "key1=value;key2=value;...".
    infos = fields[7].split(';')

    for i, info in enumerate(infos, 1):
        # info should be "key=value".
        try:
            key, value = info.split('=')
        # But sometimes it is just "value", so we'll make our own key.
        except ValueError:
            key = 'INFO{}'.format(i)
            value = info
        # Set the value to None if there is no value.
        result[key] = _get_value(value)

    return result

def lines(filename, Chr):
    """Open an optionally gzipped VCF file and generate an OrderedDict for
    each line.
    https://gist.github.com/slowkow/6215557
    """

    #TODO: see if there is a way to first map chromosomes within file and keep this data in a temp file?

    print('opening file')
    fn_open = gzip.open if filename.endswith('.gz') else open
    with fn_open(filename, 'rt') as fh, open('Exac_parse_OUT.csv', 'w') as csvout:
        a = csv.writer (csvout)
        print ('opened file', filename)
        print ('looking for chromosome ', Chr)
        items = list(range(1, 23))
        l = len(items)
        FinalResults = OrderedDict()
        #printProgress(i, l, prefix = 'Progress:', suffix = 'Complete', barLength = 50)

        for good_line in find_good_lines(fh):

            #print('Searching within chromosome ', Chr)


            query=parse(good_line)['CSQ']



            #strip first row, dont need that data
            #query = query[1]

            #search for protein coding variants for POI
            if any(ENST in s for s in query):
                #print (query)

                #combined = good_line['AC']+query

                #print (good_line)
                #print(type(good_line))
                #print (':::::')

                p_good_line = parse(good_line)
                #print(p_good_line)
                print(':::::')
                print ((p_good_line)['AC'])
                print ((p_good_line)['AF'])


                for k, e in p_good_line.items():
                    FinalResults.setdefault(k,[]).append(e)

                #a.writerows([x.split('|') for x in query] + p_good_line['AC'] + p_good_line['AF'] )
    #print(FinalResults)

    output_dict = json.loads(json.dumps(FinalResults))
    with open ('ExacDUMP.txt', 'w') as outfile:
        json.dump(output_dict, outfile)


def find_good_lines(fh):

    for line in fh:
        if line.startswith('#'):
            continue
        if line[0]<Chr:

            continue
        if line[0]==Chr:

            yield line
        if line[0]>Chr:
            break


def _get_value(value):
    """Interpret null values and return ``None``. Return a list if the value
    contains a comma.
    """
    if not value or value in ['', '.', 'NA']:
        return None
    if ',' in value:
        return value.split(',')
    return value

def filterExACoutput(file):
    print('reading file')
    df = pd.read_csv(file)
    col_names = []

    col_names = [i for i in string.printable[:len(df.columns)]]
    df.columns = [c.rstrip() for c in df.columns]
    print (col_names)
    df.columns = [col_names]
    #df[1] = [col_names[1].astype(str)]
    #df[]


    print ('done')
    df = df[df['1']=="missense_variant"]
    df = df[df['3']==GENE]
    df = df[df['6']==ENST]
    df.to_csv('Filtered_Exac_OUT.csv')

    print (df.head())


# Print iterations progress
def printProgress (iteration, total, prefix = '', suffix = '', decimals = 1, barLength = 100):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        barLength   - Optional  : character length of bar (Int)
        http://stackoverflow.com/questions/3173320/text-progress-bar-in-the-console
    """
    formatStr       = "{0:." + str(decimals) + "f}"
    percents        = formatStr.format(100 * (iteration / float(total)))
    filledLength    = int(round(barLength * iteration / float(total)))
    bar             = '█' * filledLength + '-' * (barLength - filledLength)
    sys.stdout.write('\r%s |%s| %s%s %s' % (prefix, bar, percents, '%', suffix)),
    if iteration == total:
        sys.stdout.write('\n')
    sys.stdout.flush()

def main ():

    #findPROVEANscores(ENSP)
    #formatPROVEAN(FILENAME)

    #mineExAC(ENST, Chr)

    lines('ExAC.r0.3.1.sites.vep.vcf.gz', Chr)

    #filterExACoutput('Exac_parse_OUT.csv')


    #mineMutPred(UniProt,ENST)
    #mine_dbNSFP(Chr, ENSG)
    #extract_dbNSFP(FILENAME4)



if __name__ == '__main__':
    main()




