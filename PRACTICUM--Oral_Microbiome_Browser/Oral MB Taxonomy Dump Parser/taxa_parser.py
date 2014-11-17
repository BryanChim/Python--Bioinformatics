__author__ = 'Bryan'
import argparse

parser = argparse.ArgumentParser(description='get taxa names from prof sum, map to names/nodes.dmp')
parser.add_argument('profile_summary_file', type=str, help="name of profile summary file")
parser.add_argument('output_file', type=str, help='name of output file')

args = parser.parse_args()

psum = open(args.profile_summary_file, 'r')

header = psum.readline()
taxadict = {}

for line in psum:
    line = line.strip()
    linevars = line.split('\t')


    linevars[3] = linevars[3].strip()
    if linevars[3] == 'Root':
        linevars[3] = 'root'
    taxadict[linevars[3]] = ['','','']


psum.close()

#print taxadict

namesdmp = open ('names.dmp', 'r')

for line in namesdmp:
    line = line.strip()
    linevars = line.split('\t')
    for key in taxadict:
        if linevars[2] == key:
            taxadict[key][0] = linevars[0]

#print taxadict;

namesdmp.close()

nodesdmp = open('nodes.dmp', 'r')

for line in nodesdmp:
    line = line.strip()
    linevars = line.split('\t')
    for key in taxadict:
        if linevars[0] == taxadict[key][0]:
            taxadict[key][1] = linevars[4]
            taxadict[key][2] = linevars[2]

nodesdmp.close()

out_file = open (args.output_file, 'w')

print >> out_file, 'taxa_id\tname\tlevel\tparent_taxa_id'
for key in taxadict:
    print >> out_file, taxadict[key][0], '\t', key, '\t', taxadict[key][1], '\t', taxadict[key][2]

