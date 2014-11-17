#### Parser_GenBankFile

#    PURPOSE: Parse through a GenBank (.gbk) file with a series of regular expressions
#		 and extracts gene information.
#		
#		 
#    INPUT: via command-line as argparse argument 1: 
#		'gbk_input' - GenBank input file
#
#
#    OUTPUT: via command-line as argparse argument 2:
#		'filterout' - output file containing extracted gene information, in the format:
#			<Organism>	<Chromosome>	<Contig>	<Gene_Name>	<Strand_Orientation>	<Start>	<End>
#			...
#			...
#		NOTE: each line of data is gene-centric 
#			-- there may be multiple genes found (therefore multiple lines written) for the same contig
#

import re
import argparse

# construct ArgumentParser and add arguments for file in and file out
parser = argparse.ArgumentParser(description='Extracts contig and gene information')
parser.add_argument('gbk_input', type=argparse.FileType('r'), help="hf_ref input file")
parser.add_argument('filterout', type=argparse.FileType('w'), help="hr_ref filtered output file")

# get arguments from parser
args = parser.parse_args()

# print header line to filterout
print >> args.filterout, 'Organism\tChromosome\tContig\tGene_Name\tStrand_Orientation\tStart\tEnd\n'

# initialize variables
(organism, chromosome, contig, gene_name, str_orientation, start_location, end_location) = ("","","","","","","")
seeking_gene_name = False

# for loop, iterates through hsref_input input file
for line in args.hsref_input:

    # regex to find next Contig instance. capture and save to "contig" variable
    if re.match("LOCUS", line):
        contigC = re.search('(\w\w_\d+)\s*', line)
        contig = contigC.group()
        contig = contig.strip()
        #print contig

    # regex to find organism name. capture and save to "organism" variable
    if re.search("ORGANISM\s+", line):
        organismC = re.search("ORGANISM\s+(\w+\s\w+)\s*", line)
        organism = organismC.group(1)
        #print organism

    # regex to find chromosome. capture and save to "chromosome" variable
    if re.search("\/chromosome=.*", line):
        chromosomeC = re.search("\/chromosome=\"(.+)\"", line)
        chromosome = chromosomeC.group(1)
        #print chromosome

    # regex to find gene start/end locations for positive/+ strands
    # prints out previous gene info to filterout if it exists
    # reinitializes gene_name to the null string to prevent mismatching names and sequences
    # seeking_gene_name set to True - allows for execution of IF statement below to find the corresponding gene name
    if re.search("\s+gene\s\s+<?\d+", line):
        #print >> args.filterout, line
        if contig and organism and chromosome and gene_name:
            print >> args.filterout, organism, '\t', chromosome, '\t', contig, '\t', gene_name, '\t', str_orientation, '\t', start_location, '\t', end_location
            gene_name = ""

        str_orientation = "+"
        locationC = re.search("<?(\d+)\.\.>?(\d+)", line)
        start_location = locationC.group(1)
        end_location = locationC.group(2)
        seeking_gene_name = True

        #print line, str_orientation, start_location, end_location

    # regex to find gene start/end locations for complement/- strands
    # prints out previous gene info to filterout if it exists
    # reinitializes "gene_name" to the null string to prevent mismatching names and sequences
    # "seeking_gene_name" set to True - allows for execution of IF statement below to find the corresponding gene name
    if re.search("\sgene\s\s+complement\(<?\d+\.\.>?\d+\)", line):
        #print >> args.filterout, line
        if contig and organism and chromosome and gene_name:
            print >> args.filterout, organism, '\t', chromosome, '\t', contig, '\t', gene_name, '\t', str_orientation, '\t', start_location, '\t', end_location
            gene_name = ""

        str_orientation = "-"
        locationC = re.search("complement\(<?(\d+)\.\.>?(\d+)\)", line)
        start_location = locationC.group(1)
        end_location = locationC.group(2)
        seeking_gene_name = True

        #print line, str_orientation, start_location, end_location

    # regex to find gene name
    # must have found a new gene strand (see above) - prevent mix up of gene name and strand
    # capture and save in gene_name
    # reset "seeking_gene_name" to false afterwards
    if re.search("\/gene=.*", line) and seeking_gene_name:
        gene_nameC = re.search("\/gene=\"(.*)\"", line)
        gene_name = gene_nameC.group(1)
        seeking_gene_name = False

    # regex that checks for "ORIGIN" indicating the beginning of the contig sequence itself
    # indicates that there are no more genes for the above contig -- print out the last gene's info to filterout
    # resets all variables to null string for next round of genes
    if re.match("ORIGIN", line) and gene_name:
            print >> args.filterout, organism, '\t', chromosome, '\t', contig, '\t', gene_name, '\t', str_orientation, '\t', start_location, '\t', end_location
            organism, chromosome, contig, gene_name, str_orientation, start_location, end_location = "","","","","","",""

# close files
args.hsref_input.close()
args.filterout.close()

