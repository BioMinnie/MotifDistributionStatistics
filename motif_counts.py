#######################################################################################################
#                                                                                                     #
## AUTHOR: Melinda M Ashcroft                                                                         #
## AFFILIATION: Beatson Lab | SCMB - University of Queensland                                         #
#                                                                                                     #
## PURPOSE: Obtain count of each motif per input file specified segments							                #
#                                                                       							                #
## PUBLISHED: The following script is tailored to suit the data published in Forde et al. 2015        #
## NOTE: Thank you Areej Alsheikh for your help                                                       #                                             #
#######################################################################################################

# Import modules
from collections import Counter
import csv
import sys

# Create empty lists
motif_start_pos = []
segment_start = []
segment_end = []
motif_counts = []

# Open both input files, fin1 is list of motifs and coordinates
# fin2 is genome coordinates when split into segments (ie 1kb)
# fout is  csv file of the output
with open(sys.argv[1], 'r') as fin1, \
open(sys.argv[2], 'r') as fin2:
    with open(sys.argv[3], 'w') as fout:
        
        # Specify the csv writer features and create a header row
        writer = csv.writer(fout, delimiter='\t', quotechar='', quoting=csv.QUOTE_NONE)
        headers = ('counts',)
        writer.writerow(headers)
              
        # Specify action for first input file
        motifs_file = csv.reader(fin1, delimiter='\t', quotechar='', quoting=csv.QUOTE_NONE)
        motifs_headers = fin1.readline() # reads and ignores the header line 
        for line in motifs_file:
            motif_start_pos.append(line[1]) # coordinate is in the second column (python counts from zero)
        
        # Specify action for second input file
        segments_1kb = csv.reader(fin2, delimiter='\t', quotechar='', quoting=csv.QUOTE_NONE)
        for row in segments_1kb:                # iterates through the list of 1kb segments
            segment_start.append(row[1])        # stores start position of 1kb segment
            segment_end.append(row[2])          # stores end position of 1kb segment
        
        # Specify action for output file
        for i in range(len(segment_start)):
            counts = 0
            for j in range(len(motif_start_pos)):
                if (int(motif_start_pos[j]) >= int(segment_start[i])) and (int(motif_start_pos[j]) <= int(segment_end[i])):
                    counts = counts + 1
            writer.writerow([counts])
