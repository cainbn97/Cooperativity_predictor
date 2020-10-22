'''

Reads and analyzes consensus sequence from long Homer motifs


'''
import os
import re
import operator
import numpy as np
import argparse

## parse arguments from command line
parser = argparse.ArgumentParser(description='Defines consensus sequence')
parser.add_argument('--target', help='Define TF you would like to analyze', 
    required=True)
parser.add_argument('--l_motif', type=int, help='Define TF you would like to analyze', 
    required=True)
args = parser.parse_args()

target = args.target
motif_number = args.l_motif


## read consensus sequence
c = 0
path = os.getcwd()
file = path + '/' + target + '_2_homer_denovo_long/homerResults/motif' + str(motif_number) + '.motif'
with open(file,'r') as readfile:
    for line in readfile:
        if c == 0:
            seq_temp = line.split('>')[1]
            seq = seq_temp.split('\t')[0]
        c = c + 1


## Find possible dimer sites - lots of flexibility here
## -- could also do 1 perfect, 1 very flexible site

Sites = ['TAAT','ATTA']


for s in Sites:
    print('Starting search with', s)
    Start = []
    End = []
    Group = []
    nmatches = 0
    for match in re.finditer(s,seq):
        nmatches +=1
        Start.append(match.start())
        End.append(match.end())
        Group.append(match.group())
    if nmatches == 2:
        spacer = Start[1] - End[0]
        print('Dimer found using', s)
        print('Sequence is: ', seq[Start[0]:End[1]])
        print('Site 1 is: ', Group[0])
        print('Site 2 is: ', Group[1])
        break
    elif nmatches > 2:
        print(nmatches,' spacers detected')

    elif nmatches == 1:
        print('No matches found.')
        print(Start)
        print(End)
        print(Group)
        print('Sequence is: ', seq[Start[0]:End[0]])
        
        # see if less flexible second site can be found
        Sites_flex = '[W|A|T]{4}'
        for match_2 in re.finditer('[W|A|T]{4}'):
            match_2 +=1
            Start_2.append(match_2.start())
            End_2.append(match_2.end())
            Group_2.append(match_2.group())
        if nmatches > 1:
            for n in Start_2:
                if Start[0] != Start_2[n]:
                    
                    

