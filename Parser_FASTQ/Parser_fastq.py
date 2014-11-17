#### Parser_fastq
#
#    PURPOSE: Calculates summary statistics and positional information for input FASTQ reads		
#
#		
#    INPUT: via command-line as argparse argument 1:
#		'file_in' - input FASTQ file
#
#
#    OUTPUTS: via command-line as argparse arguments 2-4:
#
#		'file_out_summary' - output file for summary statistics, in the format ~
#			<Name of fastq file>
#			<Number of reads>
#			<Average read length>
#			<Minimum read length>
#			<Maximum read length>
#			<Average read quality>
#			<Percent of reads above average>
#			<Overall GC percent>
#
#		'file_out_position_quality' - output file for quality information, in the format ~
#			<position number>	<total # reads at this position>	<average Phred Quality Score>
#
#		'file_out_position_acgt' - output file for nucleotide counts and frequencies, by position ~
#			<position number> 	<total # reads at this position>	<count and frequencies for A - C - G - T - N>
#

import re
import argparse

# create ArgumentParser
parser = argparse.ArgumentParser(description='Read in gene strand data and output +/- count and gene lengths')

# add arguments for file in and file outs
parser.add_argument('file_in', type=argparse.FileType('r'), help="Name of read-in file")
parser.add_argument('file_out_summary', type=argparse.FileType('w'), help="Name of summary out file")
parser.add_argument('file_out_position_quality', type=argparse.FileType('w'), help="Name of out position_quality out file")
parser.add_argument('file_out_position_acgt', type=argparse.FileType('w'), help="Name of position_acgt out file")

# store command line arguments into args
args = parser.parse_args()

# initialize variables
number_of_reads = 0
minimum_read_length = 1000000
maximum_read_length = 0
total_cumulative_length = 0
current_read_quality = 0
GC_count = 0

# initialize dictionaries
read_quality = {}
number_of_reads_this_position = {}
nucl_reads_by_position = {}
cumulative_quality_at_position = {}

# initialize iterators
current_position = 1
current_sequence = 1

# read in input file
for line in args.file_in:

    # check for beginning of a new read, with @ at beginning
    if re.match("^@", line):

        # increment read number
        number_of_reads += 1

        # go to next line, containing sequence, strip whitespace, store length
        line = args.file_in.readline()
        line = line.strip()
        length = len(line)

        # add to whole-file cumulative length
        total_cumulative_length += length

        # check for and update max/min read lengths
        if length > maximum_read_length:
            maximum_read_length = length
        if length < minimum_read_length:
            minimum_read_length = length

        # iterate through nucleotides in the current read
        while current_position <= length:

            # increment/initialize counter for reads at current position
            if current_position in number_of_reads_this_position:
                number_of_reads_this_position[current_position] += 1
            else:
                number_of_reads_this_position[current_position] = 1

            # initialize/increment nucleotide counts at current position
            # increment GC count if applicable
            for char in ('A', 'C', 'G', 'T', 'N'):
                if (current_position, char) not in nucl_reads_by_position:
                    nucl_reads_by_position[(current_position, char)] = 0

            if line[current_position - 1] == 'A':
                nucl_reads_by_position[(current_position, 'A')] += 1

            if line[current_position - 1] == 'C':
                GC_count += 1
                nucl_reads_by_position[(current_position, 'C')] += 1

            if line[current_position - 1] == 'G':
                GC_count += 1
                nucl_reads_by_position[(current_position, 'G')] += 1

            if line[current_position - 1] == 'T':
                nucl_reads_by_position[(current_position, 'T')] += 1

            if line[current_position - 1] == 'N':
                nucl_reads_by_position[(current_position, 'N')] += 1

            # increment current_position, reiterate through while loop
            current_position += 1

        # reset current_position to 1, for iterating through the quality line
        current_position = 1

        # initialize/reset summed read quality to 0
        current_total_read_quality = 0

        # skip to the quality line and strip its whitespace
        line = args.file_in.readline()
        line = args.file_in.readline()
        line = line.strip()

        # iterate through quality line, calculate Phred QScore for each position - update dictionaries
        while (current_position <= length):
            current_read_quality = (ord(line[current_position - 1]) - 33)
           # print ('ReadQuality:', current_read_quality)

            if current_position in cumulative_quality_at_position:
                cumulative_quality_at_position[current_position] += current_read_quality
            else:
                cumulative_quality_at_position[current_position] = current_read_quality

            current_total_read_quality += current_read_quality

            current_position += 1

        # after iterating through whole line, calculate average read quality for the current sequence
        read_quality[current_sequence] = (current_total_read_quality/length)

        # reset current_position back to 1 for the next sequence read in
        current_position = 1

        # increment counter for sequence (to be used as the next key in the read_quality dictionary)
        current_sequence += 1


# calculate average read length and GC percentage
avg_read_length = (total_cumulative_length/number_of_reads)
GC_percentage = 100*(GC_count/total_cumulative_length)

# calculate average read quality of whole data set
cumulative_read_quality = 0
for key in read_quality:
    cumulative_read_quality += read_quality[key]
average_read_quality = (cumulative_read_quality/number_of_reads)

# calculate percent of reads above average
count_for_above_average_quality = 0
for key in read_quality:
    if read_quality[key] > average_read_quality:
        count_for_above_average_quality += 1
percent_of_reads_above_average = 100*(count_for_above_average_quality/number_of_reads)

# write data to the summary file
args.file_out_summary.write('Name of fastq file: %r\n' % (args.file_in.name))
args.file_out_summary.write('Number of reads: %r\n' % (number_of_reads))
args.file_out_summary.write('Average read length: %r\n' % (avg_read_length))
args.file_out_summary.write('Minimum read length: %r\n' % (minimum_read_length))
args.file_out_summary.write('Maximum read length: %r\n' % (maximum_read_length))
args.file_out_summary.write('Average read quality: %r\n' % (average_read_quality))
args.file_out_summary.write('Percent of reads above average: %r\n' % (percent_of_reads_above_average))
args.file_out_summary.write('Overall GC percent: %r\n' % (GC_percentage))

# write data to the position_quality file
args.file_out_position_quality.write('Position Number\t# Reads at this Position\tAvg Quality at this Position\n')
for pos in range(1, maximum_read_length + 1):
    args.file_out_position_quality.write('%s\t%s\t%s\n' % (pos, number_of_reads_this_position[pos], '{0:.2f}'.format(cumulative_quality_at_position[pos]/number_of_reads_this_position[pos])))

args.file_out_position_acgt.write('Position Number\t# Reads at this Position\t#A\t%A\t#C\t%C\t#G\t%G\t#T\t%T\t#N\t%N\n')

# write data to the position_acgt file
for pos in range(1, maximum_read_length + 1):
    args.file_out_position_acgt.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' %
                                     (pos,
                                      number_of_reads_this_position[pos],
                                      nucl_reads_by_position[pos, 'A'],
                                      "{0:.2f}".format(100*nucl_reads_by_position[pos, 'A']/number_of_reads_this_position[pos]),
                                      nucl_reads_by_position[pos, 'C'],
                                      "{0:.2f}".format(100*nucl_reads_by_position[pos, 'C']/number_of_reads_this_position[pos]),
                                      nucl_reads_by_position[pos, 'G'],
                                      "{0:.2f}".format(100*nucl_reads_by_position[pos, 'G']/number_of_reads_this_position[pos]),
                                      nucl_reads_by_position[pos, 'T'],
                                      "{0:.2f}".format(100*nucl_reads_by_position[pos, 'T']/number_of_reads_this_position[pos]),
                                      nucl_reads_by_position[pos, 'N'],
                                      "{0:.2f}".format(100*nucl_reads_by_position[pos, 'N']/number_of_reads_this_position[pos])))


# end of Parser_fastq.py