__author__ = 'Bryan'
import collections
vardict = {}

cdf = open ('OralMicrobiome_Clinical_Mapping_Data.txt', 'r')
cdout = open ('SAMPLEATTRIBUTES_loadData.dat','w')
header = cdf.readline()

#print header

variables = header.split('\t')
variables.pop(0)
for var in variables:
    var = var.strip()
variables[33] = 'PI-S'
#print len(variables)
#print variables;

for line in cdf:
    line = line.strip()
    linevars = line.split('\t')
    vardict[linevars[0]] = linevars

    #print vardict[linevars[0]]
    linevars.pop(0)


ovardict = collections.OrderedDict(sorted(vardict.items()))
cdf.close()
#print len(vardict['1_D_A01_2010_12_17_2'])

print >> cdout, 'sample','\t','attribute','\t','value'

for key in ovardict:
    for i in range(len(variables)):
        print >> cdout, key,'\t',variables[i],'\t',ovardict[key][i]


cdout.close()