'''
Missense variant miner:
find missense variants in ExAC, score them using PROVEAN, MutPred and dbNSFP
return csv file.



'''
from __future__ import print_function
#from chardet.universaldetector import UniversalDetector
import sys, gzip, tarfile, os
import pandas as pd
import csv
import zipfile
from collections import OrderedDict
import string
import json

VCF_HEADER = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']

# Sets protein ID to search in dataframe

#ENST = "ENST00000448921" #A1AT
#ENSP = "ENSP00000348068" #A1AT
#ENSG = "ENSG00000197249" #A1AT

#ENSP = "ENSP00000003084" #CFTR
#ENSP = "ENSP00000269228" #NPC1
#ENSP = "ENSP00000262304" #PKD1
ENSP = "ENSP00000262410" #MAPT
#ENSG = "ENSG00000001626" #CFTR
#ENSG = "ENSG00000141458" #NPC1
#ENSG = "ENSG00000008710" #PKD1
#ENSG = "ENSG00000118762" #PKD2
#ENSG = "" #PKD2
#ENSG = "ENSG00000186868" #ENSG
ENSG = "ENSG00000186868" #MAPT
# ENSG = "ENSG00000272636" #Diagnostic - beginning of Chr17
#ENST = "ENST00000003084" CFTR
#ENST = 'ENST00000269228' #NPC1i
#ENST = "ENST00000262304" #PKD1
#ENST = "ENST00000237596" #PKD2
ENST = "ENST00000262410" #MAPT


GENE = "MAPT"
FILENAME = "CFTR_PROV_extract.csv"
FILENAME1 = "CFTR_PROVEANScores.csv"
FILENAME2 = "A1AT_ExACScores.csv"
FILENAME3 = "A1AT_MutPredScores.csv"
FILENAME4 = "MAPT_dbNSFPa_output.csv"
FILENAME5 = "dbNSFP_extract.csv"
FILENAME6 = "A1AT_ExACScores.csv"
#UniProt = "P13569" CFTR
#UniProt = "P01009" #A1AT
#UniProt = "O15118" #NPC1
#UniProt = "P98161" #PKD1
#UniProt = "P10636" #MAPT
#Chr = '18' NPC1
#Chr = "14" #A1AT
#Chr = "4" #PKD2
#Chr = "16" #PKD1
#Chr = "1"

#FILENAME6 = "SOD1_ExACScores.csv"
#ENSG = "ENSG00000142168" #SOD1
#ENST = "ENST00000270142" #SOD1
#Chr = "21" #SOD1
Chr = "17" #MAPT


# change directory to working with DAta
os.chdir("../Data/")
cwd = os.getcwd()


# df = pd.read_csv('foo.csv', index_col=0)
def findPROVEANscores(ENSP):
    # read from tsv.gz file
    with gzip.open('PROVEAN_scores_ensembl66_human.tsv.gz', 'rt') as tsvin, open(FILENAME1, 'wt') as csvout:
        csvout = csv.writer(csvout)
        tsvin = csv.reader(tsvin, delimiter='\t', quoting=csv.QUOTE_NONE)
        for i in range(1):
            row1 = next(tsvin)
            print(row1)
            csvout.writerows([row1])
            i = +1
        for row in tsvin:
            count = row[0]
            # print (count, ENSP)
            if count == ENSP:
                csvout.writerows([row[0:23]])


def formatPROVEAN(input):
    df = (pd.read_csv(input))
    short_df = df.drop(df.columns[[0, 1]], axis=1)
    short_df['AA'] = (short_df.T.idxmax())
    short_df['position'] = df['position'].apply(str)
    short_df['Mutation'] = short_df['AA'].astype(str) + short_df['position'].astype(str)
    df.to_csv(FILENAME, sep='\t')
    print(short_df)


def dataframeExAC(SYMBOL):
    print("Open")

    # read from tsv.gz file
    with gzip.open('ExAC.r0.3.1.sites.vep.vcf.gz', 'rt') as tsvin, open(FILENAME2, 'wt') as csvout:
        print("open")
        comments = count_comments(tsvin)
        result = OrderedDict()
        csvout = csv.writer(csvout)
        tsvin = pd.read_table(tsvin, skiprows=comments, names=VCF_HEADER,
                              usecols=range(8))  # delimiter='|',quoting=csv.QUOTE_NONE, skiprows = comments)
        print('looking for this query: ', SYMBOL)

        for i, line in enumerate(lines(tsvin)):
            for key in line.keys():
                if key not in result:
                    result[key] = [None] * i

            for key in result.keys():
                result[key].append(line.get(key, None))

        return pd.DataFrame(result)


def generatorExAC(filename, Chr):
    print('searching Chrom:', Chr, 'in file: ', filename)
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

def mineMutPred(ENST, UniProt, Chr):
    ChrMutPredFilesDict = {
        '1': 'dbnsfp_chr1_mutpred.txt',
        '2': 'dbnsfp_chr2_mutpred.txt',
        '3': 'dbnsfp_chr3_mutpred.txt',
        '4': 'dbnsfp_chr4_mutpred.txt',  # ENSG is row 19"
        '5': 'dbnsfp_chr5_mutpred.txt',
        '6': 'dbnsfp_chr6_mutpred.txt',
        '7': 'dbnsfp_chr7_mutpred.txt',
        '8': 'dbnsfp_chr8_mutpred.txt',
        '9': 'dbnsfp_chr9_mutpred.txt',
        '10': 'dbnsfp_chr10_mutpred.txt',
        '11': 'dbnsfp_chr11_mutpred.txt',
        '12': 'dbnsfp_chr12_mutpred.txt',
        '13': 'dbnsfp_chr13_mutpred.txt',
        '14': 'dbnsfp_chr14_mutpred.txt',
        '15': 'dbnsfp_chr15_mutpred.txt',
        '16': 'dbnsfp_chr16_mutpred.txt',
        '17': 'dbnsfp_chr17_mutpred.txt',
        '18': 'dbnsfp_chr18_mutpred.txt',
        '19': 'dbnsfp_chr19_mutpred.txt',
        '20': 'dbnsfp_chr20_mutpred.txt',
        '21': 'dbnsfp_chr21_mutpred.txt',
        'X': 'dbnsfp_chrX_mutpred.txt',
        'Y': 'dbnsfp_chrY_mutpred.txt',
    }

    # TODO: This is working but is very slow, need to make into Generator?
    # read from tsv.gz file
    with tarfile.open('Mutpred.tar.gz') as tsvin, open(FILENAME3, 'wt') as csvout:
        a = csv.writer(csvout)
        first_row = ('Genomic_location', 'ProteinID', 'MUTATION', 'MUTPRED Score', 'Hypothesis')
        a.writerows([first_row])
        tar = tsvin.extractfile(ChrMutPredFilesDict[Chr])
        print(tar)

        tp = pd.read_table(tar, delimiter='\t', quoting=csv.QUOTE_NONE)
        print(type(tp))

        # tsvin = csv.reader(tsvin, delimiter='\t',quoting=csv.QUOTE_NONE)
        print('looking for this query: ', ENST, UniProt, 'in chromosome',Chr)
        # TODO: write first row to file
        # for i in range(1):
        #     row1 = next(tp)
        #     print(row1)
        #     print('found ', len(row1), 'rows' )
        #     csvout.writerows([row1])
        #     i=+1
        variants = 0
        for index, row in tp.iterrows():

            count = row[1]
            # print(count,ENST, UniProt)
            if count == ENST or count == UniProt:
                variants += 1
                print(variants, ' writing variant scores ', row[0:2])
                a.writerows([row[0:len(row)]])

        print('found ', variants, 'variant MutPred Scores')


def mine_dbNSFP(Chr, ENSG):
    # TODO: should be able to read file names and iterate through them
    chrfilesdict = {
        '1': 'dbNSFP3.2a_variant.chr1',
        '2': 'dbNSFP3.2a_variant.chr2',
        '3': 'dbNSFP3.2a_variant.chr3',
        '4': 'dbNSFP3.2a_variant.chr4',
        '5': 'dbNSFP3.2a_variant.chr5',
        '6': 'dbNSFP3.2a_variant.chr6',
        '7': 'dbNSFP3.2a_variant.chr7',
        '8': 'dbNSFP3.2a_variant.chr8',
        '9': 'dbNSFP3.2a_variant.chr9',
        '10': 'dbNSFP3.2a_variant.chr10',
        '11': 'dbNSFP3.2a_variant.chr11',
        '12': 'dbNSFP3.2a_variant.chr12',
        '13': 'dbNSFP3.2a_variant.chr13',
        '14': 'dbNSFP3.2a_variant.chr14',
        '15': 'dbNSFP3.2a_variant.chr15',
        '16': 'dbNSFP3.2a_variant.chr16',
        '17': 'dbNSFP3.2a_variant.chr17',
        '18': 'dbNSFP3.2a_variant.chr18',
        '19': 'dbNSFP3.2a_variant.chr19',
        '20': 'dbNSFP3.2a_variant.chr20',
        '21': 'dbNSFP3.2a_variant.chr21',
        'X': 'dbNSFP3.2a_variant.chrX',
        'M': 'dbNSFP3.2a_variant.chrM',
    }
    # read from tsv.gz file
    with zipfile.ZipFile('dbNSFPv3.2a.zip', 'r') as tsvin, open(FILENAME4, 'wt') as csvout, open('tmp.csv', 'wt') as csvout_temp:
        tp = pd.read_csv(tsvin.open(chrfilesdict['1']), delimiter='\t', quoting=csv.QUOTE_NONE, iterator=True,
                         dtype=object, chunksize=1)  # header = None)

        writer = csv.writer(csvout, dialect='excel')
        writer_temp = csv.writer(csvout_temp, dialect = 'excel')
        for i in range(1):
            row1 = next(tp)
            writer.writerows([row1])
            i+=1
        print(chrfilesdict[Chr])
        print(ENSG)
        for line in enumerate(tsvin.open(chrfilesdict[Chr])):
            line = line[1].decode("utf-8").rstrip().split('\t')
            if ENST == line[20]:
                writer.writerows([line])
            elif ENST in line[20] and ';' in line[20]:

                for i, val in enumerate(line):
                    if ';' in val:

                        val = val.split(';', 1)[0]
                        line[i] = val
                    if val == '.':


                        line [i] = '0'

                writer.writerows([line])

def db_NSFP_iterate(fh):
    #TODO: rewrite to account for Chr vs int variables and to reduce unnecessary searching.
    i='0'
    found = 0

    for line in fh:
        print (line)

        ENSG_match = line[20]

        if ENSG_match == ENSG:
            print (line)
            yield line
            found = 1
        else:
            if ENSG_match != ENSG and found ==1:
                print ("no more to be found")
                break
            continue

def cleanup_dbNSFP_extract(file):

    with open(file, 'rt') as tsvin, open(FILENAME5, 'wt') as csvout:
        dict = []
        df = pd.read_csv(tsvin, delimiter=',', encoding="utf-8-sig")

        for row in enumerate(df['FATHMM_score']):

            FAS = row[1]

            print (FAS)
            if FAS[0] == FAS[1]:
                id = row[0]
                print(df.iloc[[id]])


def extract_dbNSFP(file):

    """ There are rows with multiple comma separated values in them
        this function works to convert such rows into rows which have
        one value per line
        Problems is that it will not run unless rows have exact number of comma sep. values
        so next step is to subset these rows and then merge with the rest of the dataframe
    """
    with open(file, 'rt') as tsvin, open(FILENAME5, 'wt') as csvout:
        dict = []
        df = pd.read_csv(tsvin, delimiter=',', encoding="utf-8-sig")

        #print (df.head(1))
        #print (df['Ensembl_transcriptid'])
        #column_names = ['Ensembl_transcriptid', 'Ensembl_proteinid',  'MutationTaster_score','MutationTaster_pred', 'MutationTaster_AAE', 'FATHMM_score', 'FATHMM_pred']
        #print(column_names)
        df = df[['Ensembl_transcriptid', 'Ensembl_proteinid',  'MutationTaster_score','MutationTaster_pred', 'MutationTaster_AAE', 'FATHMM_score', 'FATHMM_pred']]

        for col in ['Ensembl_transcriptid', 'Ensembl_proteinid',  'MutationTaster_score','MutationTaster_pred', 'MutationTaster_AAE', 'FATHMM_score', 'FATHMM_pred']:
            df [col] = df[col].str.split(';')



            #print (len(row['Ensembl_proteinid']))'Ensembl_proteinid'

'Ensembl_proteinid'
        #
        # #print (df['Ensembl_transcriptid'])
        # #df.to_csv(csvout)
        # i =  df['Ensembl_transcriptid'].map(len)
        # j = np.repeat(np.arange(len(df)),i)
        # k = np.concatenate(list(map(np.arange, i)))
        # df = df.iloc[j]
        # print (df['Ensembl_transcriptid'])
        # for col in ['Ensembl_transcriptid', 'Ensembl_proteinid',  'MutationTaster_score','MutationTaster_pred', 'MutationTaster_AAE', 'FATHMM_score', 'FATHMM_pred']:
        #     df [col] = list(map(lambda xs, i: xs[i], df[col], k))
        # df.to_csv(csvout)


def count_comments(filename):
    """Count comment lines (those that start with "#") in an optionally
	gzipped file.
	:param filename:  An optionally gzipped file.
	https://gist.github.com/slowkow/6215557
	"""
    comments = 0
    print("comm")
    # fn_open = gzip.open if filename.endswith('.gz') else open
    # with fn_open(filename) as fh:
    for line in filename:
        if line.startswith('#'):
            comments += 1
            print(line)
        else:
            print("no more")
            break
    print(comments)
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
    # TODO: see if there is a way to first map chromosomes within file and keep this data in a temp file?
    print('opening file')
    fn_open = gzip.open if filename.endswith('.gz') else open
    with fn_open(filename, 'rt') as fh, open(FILENAME6, 'w') as csvout:
        a = csv.writer(csvout)
        first_row = ('GENE_ID', 'TRANSCRIPT', 'TRANSCRIPT CHANGE', 'PROTEIN CHANGE', 'AA_POS', 'AA_CHANGE', 'MUTATION',
                     'ALLELE COUNT', 'ALLELE FREQUENCY')
        a.writerows([first_row])
        print('opened file', filename)
        print('looking for chromosome ', Chr)
        items = list(range(1, 23))
        l = len(items)
        FinalResults = OrderedDict()
        # printProgress(i, l, prefix = 'Progress:', suffix = 'Complete', barLength = 50)
        for good_line in find_good_lines(fh):
            # parse data line to find protein ID
            # dict item ['CQS'] contains ENSP IDs and mutation ID in tricode
            query = parse(good_line)['CSQ']

            # search for protein coding variants for POI
            if any(ENST in s for s in query):
                # get the rest of the data for matching lines
                p_good_line = parse(good_line)
                for k, e in enumerate(query):
                    for line in e.splitlines():
                        if ENSP in line and "missense_variant" in line:  # TODO: add deletions
                            AlleleCount = ()
                            AlleleFrequency = ()
                            print(':::::')
                            # print(p_good_line['AF'])
                            print(len(p_good_line['AC']))
                            if len(p_good_line['AC']) > 1:
                                print('AC1', (p_good_line['AC']))
                                print('AF1', (p_good_line['AF']))
                                try:
                                    AlleleCount = sum(list(map(int, p_good_line[
                                        'AC'])))  # list added to break object, if iterate dont need list(
                                    AlleleFrequency = sum(list(map(float, p_good_line['AF'])))
                                except ValueError:
                                    print("error")
                                    AlleleCount = p_good_line['AC']
                                    AlleleFrequency = float(p_good_line['AF'])
                                    print('AF2:', p_good_line['AF'])
                                    print('AC2:', p_good_line['AC'])
                            else:
                                AlleleCount = int(p_good_line['AC'])
                                AlleleFrequency = float(p_good_line['AF'])

                            match = line.split( '|')  # GeneSym:match[3], geneID:match[4], ENST:match[6], ENST+change: match[10], ESNP_change:match [11], AApos:match[14], AA_change: match[15]
                            # print(match[3], match[6],match[10],match[11],match[14], match[15])
                            aa_change = match[15].split('/')
                            variant = [aa_change[0], match[14], aa_change[1]]
                            print(variant)
                            # print(':::::::::::::::::::::::::::::::::::::::::::::::')
                            row_line = [match[3], match[6], match[10], match[11], match[14], match[15],
                                        ''.join(variant), AlleleCount, AlleleFrequency]
                            print(row_line)
                            a.writerows([row_line])
    # print(FinalResults)

    output_dict = json.loads(json.dumps(FinalResults))
    with open('ExacDUMP.txt', 'w') as outfile:
        json.dump(output_dict, outfile)


def find_good_lines(fh):
    #TODO: rewrite to account for Chr vs int variables and to reduce unnecessary searching.
    i='0'
    for line in fh:

        chrom = line[0:2]
        #print (chrom, Chr)
        if line.startswith('#'):
            continue
        if chrom == Chr:


            yield line
        else:
            if line[0:2] != i:

                print (line[0:2])
                i = line[0:2]


            continue


        # if chrom > Chr:
        #     print ('end')
        #     break


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
    print(col_names)
    df.columns = [col_names]
    # df[1] = [col_names[1].astype(str)]
    # df[]


    print('done')
    df = df[df['1'] == "missense_variant"]
    df = df[df['3'] == GENE]
    df = df[df['6'] == ENST]
    df.to_csv('Filtered_Exac_OUT.csv')

    print(df.head())




def main():

    # findPROVEANscores(ENSP)
    # formatPROVEAN(FILENAME)

    #mineExAC(ENST, Chr)

    #lines('ExAC.r0.3.1.sites.vep.vcf.gz', Chr)


    #filterExACoutput('Exac_parse_OUT.csv')


    #mineMutPred(UniProt,ENST,Chr)
    mine_dbNSFP(Chr, ENSG)
    #extract_dbNSFP(FILENAME4)
    #cleanup_dbNSFP_extract(FILENAME4)


if __name__ == '__main__':
    main()
