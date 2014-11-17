#### Parser_BLAST_refseq
#
#		*** SEE ACCOMPANYING README.TXT FILE FOR DETAILED INFORMATION ***
#
#    PURPOSE: 
#    INPUT:
#    OUTPUT: 
#			


import re
import argparse
import os


# created argument parser and add arguments - see help= for descriptions
parser = argparse.ArgumentParser(description='Read in blast output file, write out filtered results according to threshold.')

parser.add_argument('refseqs', type=str, help="name of query refseq id file")
parser.add_argument('qfna', type=str, help="name of query fasta file")
parser.add_argument('dbfna', type=str, help="name of database fasta file")
parser.add_argument('blastout', type=str, help="desired name of blast output file")
parser.add_argument('filterout', type=str, help="desired name of file for output filtered by threshold")

# method for a restricted float type, for use in threshold argument
def restricted_float(t):
    t = float(t)
    if (t < 0.0 or t > 1.0):
        raise argparse.ArgumentTypeError("threshold value must be between 0.0 and 1.0" % t)
    return t

parser.add_argument('threshold', type=restricted_float, help="threshold for hit extraction")


# parse arguments into args
args = parser.parse_args()

# open refseq id file, query fasta file and file to write extracted sequences to
refseq_ids = open(args.refseqs, 'r')
query_fna = open(args.qfna, 'r')
refseq_extr = open("extracted_qrefseqs.fna", 'w')

# initialize list in which refseq ids are to be added
refseq_list = []

# iterate through refseq id file, get ids and append to refseq_list
for line in refseq_ids:
    refseq = line.strip()
    refseq_list.append(refseq)

#print refseq_list

seq_found = False

# parse through query fasta file, extract those sequences that specified in the refseq list
for line in query_fna:
        if not re.match(">gi.*", line) and seq_found:
            line = line.strip()
            print >> refseq_extr, line
        else:
            seq_found = False
            for id in refseq_list:
                if re.search(id, line):
                    line = line.strip()
                    print >> refseq_extr, line
                    seq_found = True

# close files
refseq_ids.close()
refseq_extr.close()
query_fna.close()

# blast commands - create database and query the extracted query sequences against it
# -outfmt 6: customized tab-delimited output, displays:
# query id, subject id, query length, subject length, alignment length, bitscore and evalue
os.system("makeblastdb -in"
          " " + args.dbfna + " -dbtype nucl -out blast_db")
os.system("blastn -db blast_db -query extracted_qrefseqs.fna -outfmt \"6 qseqid sseqid qlen slen length bitscore evalue\" -out " + args.blastout + " -evalue 1e-05")

# open blast output file for reading
# open filtered output file for writing
blast_output = open(args.blastout, "r")
filter_out = open(args.filterout, 'w')

# write the header to the filtered output file
filter_out.write('Mouse Refseq-id\tHuman Refseq-id\tLength of Mouse Gene\tLength of Human Gene\tAlignment-length\tBit-score\tE-score\n')

# iterate through the hits in the blast_output
for line in blast_output:

    # strip line and place contents in array
    line = line.strip()
    linearray = line.split()

    # calculates alignment length over query length, and compares it against the inputted threshold
    # if the hit meets the threshold, extract ref seq ids and other hit info -- write it all to the filtered output file
    if float(linearray[4])/float(linearray[2]) > args.threshold:

        mouse_refseqIDC = re.search("(\w\w_\d+\.\d)", linearray[0])
        mouse_refseqID = mouse_refseqIDC.group()

        human_refseqIDC = re.search("(\w\w_\d+\.\d)", linearray[1])
        human_refseqID = human_refseqIDC.group()

        filter_out.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' %
            (mouse_refseqID, human_refseqID, linearray[2],
            linearray[3], linearray[4], linearray[5],
            linearray[6]))

# close blast file and output file
blast_output.close()
filter_out.close()

# regexes designed for the default blast output
"""
for line in args.blastout:
    if re.match("Query= ", line):
        if (mrefseqID and qlength and hlength and hrefseqID):
            print mrefseqID,qlength,hrefseqID,hlength
            mrefseqID,qlength,hrefseqID,hlength = "","","",""
        mrefseqC = re.search("(NM_\d+\.\d)", line)
        mrefseqID = mrefseqC.group()
        seek_query_length = True

    if re.match("Length=\d+", line):
        lengthC = re.match("Length=(\d+)", line)
        if seek_query_length:
            qlength = lengthC.group()
            seek_query_length = False
        else:
            hlength = lengthC.group()

    if re.match("> gi|.*", line):
        if mrefseqID and qlength and hlength and hrefseqID:
            print mrefseqID,qlength,hrefseqID,hlength
            hrefseqID, hlength = "",""

        hrefseqC = re.search("(\w\w_\d+\.\d)", line)
        if hrefseqC:
            hrefseqID = hrefseqC.group()

    if re.match("Gap Penalties", line):
        if mrefseqID and qlength and hlength and hrefseqID:
            print mrefseqID,qlength,hrefseqID,hlength
"""


