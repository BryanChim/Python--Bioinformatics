__author__ = 'Bryan'

import re

threshold = 0.6
method_id = 1
total_reads = 0
previous_score_index = 0

clin = open('S1_1_2011_04_27_1_CACTGA_classified.txt', 'r')
print 'SampleID\tMethod-id\tReadID\tTaxa-Name\tTaxa-Level\tScore\n'
taxa_dict = {}
taxa_name = ''

for line in clin:

    line = line.strip()
    line = line.replace('"', '')
    linearray = line.split('\t')

    read_sample_idC = re.match("(.*?)\|\d+\|(.*?)\|.*", linearray[0])
    sample_id = read_sample_idC.group(1)
    read_id = read_sample_idC.group(2)

    current_index = linearray.index('Root')
    current_index = current_index + 2
    #previous_score_index = current_index - 3

    while current_index < len(linearray):
        if (re.match("\d\.\d+", linearray[current_index])):
            if float(linearray[current_index]) >= threshold:

               # for ind in range(previous_score_index + 1, current_index - 1):
                #    taxa_name = taxa_name + linearray[ind]
                taxa_name = linearray[current_index-2]
                print sample_id, '\t', method_id, '\t', read_id, '\t', taxa_name, '\t', linearray[current_index - 1], '\t', linearray[current_index], '\n'

                if taxa_name in taxa_dict:
                    taxa_dict[taxa_name][4] += 1
                    taxa_dict[taxa_name][5] += float(linearray[current_index])
                else:
                    taxa_dict[taxa_name] = [sample_id, method_id, taxa_name, linearray[current_index - 1], method_id, float(linearray[current_index])]

                #previous_score_index = current_index
                current_index = current_index + 3
                taxa_name = ''

            else:
                break
        else:
            current_index += 1


    total_reads = total_reads + 1

clin.close()

for taxa in taxa_dict:
    percent_of_total = float(taxa_dict[taxa][4]) / float(total_reads) * 100
    avg_score = taxa_dict[taxa][5] / float(taxa_dict[taxa][4])
    print taxa_dict[taxa][0], '\t', taxa_dict[taxa][1], '\t', taxa_dict[taxa][2], '\t', taxa_dict[taxa][3], '\t', \
        taxa_dict[taxa][4], '\t', "%.2f" % percent_of_total, '\t', "%.2f" % avg_score
