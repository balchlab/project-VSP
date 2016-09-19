'''
Converts ExAC missense protein notation to single letter notation
TODO: get PROVEAN.

input: .txt file with "Arg13Pro" style notation
output: .txt file with  "R13P" notation
'''

import sys, traceback,subprocess,gzip,glob, tarfile,os, signal
import pandas as pd

import s3fs

fs = s3fs.S3FileSystem()

with fs.open('s3://bucket_name/objkey') as f:
    df = pd.read_csv(f, compression='gzip', nrows=5)



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

#
# def Get_PROVEAN_Score():
#     os.chdir("../Data")
#     tar = tarfile.open("test.tar")
#     for member in tar.getmembers():
#         f=tar.extractfile(member)
#         content=f.read()
#         print ("%s has %d newlines" %(member, content.count("\n")))
#         print ("%s has %d spaces" % (member,content.count(" ")))
#         print ("%s has %d characters" % (member, len(content)))
#     sys.exit()
# tar.close()

def filter_from_tsvs(input_dir,string):

    tsvs = glob.glob(os.path.join(input_dir,'*.tsv*'))
    open_cmd=open
    for tsvfile in tsvs:
        print (os.path.splitext)
        extension = os.path.splitext(tsvfile)[1]
        if extension == ".gz":
          open_cmd = gzip.open
    print (open_cmd)
    try:
        print (subprocess.check_output('grep string tsvfile', shell=True))
        output = subprocess.check_output('grep string tsvfile', shell=True, preexec_fn=lambda: signal.signal(signal.SIGPIPE, signal.SIG_DFL))
        print(output)

    except Exception as e:
        print ("%s" %e)
        print ("%s" %traceback.format_exc())




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