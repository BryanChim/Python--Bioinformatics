__author__ = 'Bryan'
import re
import argparse
import os

read_dict = {}


parser = argparse.ArgumentParser(description='Get reads and read info')
parser.add_argument('dir_of_fasta_files', type=str, help="directory to input fasta files")
parser.add_argument('output_file', type=argparse.FileType('w'), help="output file of read data")

args = parser.parse_args()

fasta_files = os.listdir(args.dir_of_fasta_files)

for file in fasta_files:


#parser.add_argument('dir_of_fasta_files', type=str, help="directory to input fasta file")

#path_to_current_fasta = ''
# get arguments from parser

    fa_in = open(file, 'r')
    sample_id_from_fileC = re.search("(.*?)_[ATCG]+?_(?:read|test).*", file)
    sample = sample_id_from_fileC.group(1)

    #fasta_files = os.listdir(args.dir_of_fasta_files)

    for line in fa_in:
        line = line.strip()
        if (re.match('>.*', line)):
            read_infoG = re.match('>.*\|\d+\|(.*?)\|(\d+)\|(\d+)', line)
            name = read_infoG.group(1)
            length = read_infoG.group(2)
            qualityscore = read_infoG.group(3)

            print >> args.output_file, name, '\t', sample, '\t', length, '\t', qualityscore

    fa_in.close()

args.output_file.close()

        #read_dict[name] = (length, qualityscore)



