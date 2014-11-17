__author__ = 'Bryan'
import re
import argparse
import os

#construct ArgumentParser and add arguments for file in and file out
parser = argparse.ArgumentParser(description='Extracts taxa classification information')
parser.add_argument('path_of_classifierjar', type=str, help="path to RDP classifier.jar file")
parser.add_argument('dir_of_fasta_files', type=str, help="directory to input fasta file")

path_to_current_fasta = ''
# get arguments from parser
args = parser.parse_args()

fasta_files = os.listdir(args.dir_of_fasta_files)

for file in fasta_files:

    sample_id_from_fileC = re.search("(.*?)_[ATCG]+?_(?:read|test).*", file)
    sample_id_from_file = sample_id_from_fileC.group(1)
    path_to_current_fasta = args.dir_of_fasta_files + '\\' + file

    #print sample_id_from_file

    os.system("java -Xmx1g -jar \"" + args.path_of_classifierjar + "\" classify -c 0.5 -o " + sample_id_from_file + "_classified.txt -h " + sample_id_from_file + "_hierarchy.txt \"" + path_to_current_fasta + "\"")

    classified_filename = sample_id_from_file + "_classified.txt"
    read_assignment_output_filename = sample_id_from_file + "_read_assignment.txt"
    profile_summary_output_filename = sample_id_from_file + "_profile_summary.txt"


    clin = open(classified_filename, 'r')
    ra_out = open(read_assignment_output_filename, 'w')
    ps_out = open(profile_summary_output_filename, 'w')

    method_id = '2'
    read_id,taxa_name,score,percent_of_total,avg_score = ('','','','','')
    threshold = 0.8
    taxa_dict = {}
    total_reads = 0
    current_index = 0
    #previous_score_index = 0

    print >> ra_out, 'SampleID\tMethod-id\tReadID\tTaxa-Name\tTaxa-Level\tScore'

    for line in clin:

        line = line.strip()
        line = line.replace('"', '')
        linearray = line.split('\t')

        read_sample_idC = re.match(".*?\|\d+\|(.*?)\|.*", linearray[0])
        read_id = read_sample_idC.group(1)

        current_index = linearray.index('Root')
        current_index = current_index + 2
        #previous_score_index = current_index - 3


        while current_index < len(linearray):
            if (re.match("\d\.\d+", linearray[current_index])):
                if float(linearray[current_index]) >= threshold:

                    #for ind in range(previous_score_index + 1, current_index - 1):
                    #    taxa_name = taxa_name + linearray[ind]
                    taxa_name = linearray[current_index-2]

                    print >> ra_out, sample_id_from_file, '\t', method_id, '\t', read_id, '\t', taxa_name,'\t',linearray[current_index-1],'\t', linearray[current_index]

                    if taxa_name in taxa_dict:
                        taxa_dict[taxa_name][4] += 1
                        taxa_dict[taxa_name][5] += float(linearray[current_index])
                    else:
                        taxa_dict[linearray[current_index-2]] = [sample_id_from_file,method_id,taxa_name,linearray[current_index-1],1,float(linearray[current_index])]

                    #previous_score_index = current_index
                    current_index = current_index + 3
                    taxa_name = ''
                else:
                    break
            else:
                current_index = current_index + 1

        total_reads = total_reads + 1

    clin.close()


    print >> ps_out, 'SampleID\tMethod-id\tTaxa-Name\tTaxa-Level\t#_of_Reads\t%_of_Total\tAvg_Score'

    for taxa in taxa_dict:
        percent_of_total = float(taxa_dict[taxa][4])/float(total_reads)*100
        avg_score = taxa_dict[taxa][5]/float(taxa_dict[taxa][4])
        print >> ps_out, sample_id_from_file,'\t',taxa_dict[taxa][1],'\t',taxa_dict[taxa][2],'\t',taxa_dict[taxa][3],'\t',taxa_dict[taxa][4],'\t',"%.2f" % percent_of_total,'\t',"%.2f" % avg_score

    clin.close()
    ra_out.close()
    ps_out.close()