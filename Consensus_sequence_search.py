'''

Reads and analyzes consensus sequence from long Homer motifs


'''
import os
import re
import operator
import numpy as np
import replace
import argparse

## parse arguments from command line
parser = argparse.ArgumentParser(description='Defines consensus sequence')
parser.add_argument('TF', type=string, help='Define TF you would like to analyze', 
    required=True)
parser.add_argument('Cycle', type=int, help='Define cycle you would  like to analyze', 
    required=True)
args = parser.parse_args()

target=args.TF
cycle=args.Cycle


## read consensus sequence
path = os.getcwd()
file = path + '/' + target + '_homer_long/homerResults/motif1.motif'
with open(file) as readline:
    for line in readfile:
        if line == 0:
            seq_temp = line.split('>')[1]
            seq = seq_temp.split('\t')[0]



#Run code for every spacer possibility
Found = False
for monomer in ['TAAT','WAAT', Found == False:

    ## find perfect monomer locations
    Start = []
    End = []
    Group = []
    for match in re.finditer('TAAT',seq):
        nmatches +=1
        Start.append(match.start)
        End.append(match.end)
        Group.append(match.group)
    if nmatches == 2:
        spacer = Start[1] - End[0]
        Found = True
    elif nmatches > 2:
        print(nmatches,' spacers detected')
    elif nmatches == 1:
        print('No matches found with' 
        



    # #Forward sequence search
    # motif_for = 'TAAT[\w]{'+str(spacer)+'}TAAT'
    # result = re.finditer(motif_for, seq)

    # #set variables
    # nmatches = 0
    # checklist_for = ''

    # #search forward reads
    # for match in result:
        # nmatches += 1
        # #Save all results to later check for duplicates
        # #'+' added to prevent searches from spanning across separate reads
        # checklist_for = checklist_for+file[(match.start()-2):(match.end()+2)]+'+'

    # #reverse sequence search
    # motif_rev = 'ATTA[G|A|T|C]{'+str(spacer)+'}ATTA'
    # result = re.finditer(motif_rev, seq)

    # #search reverse reads
    # for match in result:
        # nmatches += 1
    

    # #Check for reverse reads in forward search list
    # check_for = re.finditer(motif_rev, checklist_for)

    # #set variables
    # checkmatch = 0

    # #search for reverse motif in forward reads
    # for match in check_for:
        # checkmatch +=1
        # #print(checklist_rev[(match.start()-2):(match.end()+2)])
    
    # #save replicates
    # checkmatch_tot.append(checkmatch)

    # #subtract replicates
    # nmatches = nmatches - checkmatch

    # #add data to vector  
    # percent_dimer = nmatches/nhits * 100
    # percent_dimer = round(percent_dimer, 2)
    # final_data.append(percent_dimer)
    # spacer = spacer - 1
# print(final_data)   
