#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Detects top dimer site from Homer's de novo motif analysis
'''
import os
import numpy as np
import glob
import pandas as pd
import re


## Grab motif file
c = 0
motif_final = []
path = os.getcwd()
motif_file = path + '/*_long.motif'

## Read in motif
for file in glob.glob(motif_file):
    with open(file,'r') as readfile:
        for line in readfile:
            if c == 0:
                seq_temp = line.split('>')[1]
                seq = seq_temp.split('\t')[0]
            else:
                motif_row_final = []
                motif_row = line.replace('\n','')
                for m in np.arange(0,4):
                    motif_row_final.append(float(motif_row.split('\t')[m]))
                motif_final.append(motif_row_final)          
            c = c + 1
            
motif = pd.DataFrame(motif_final, columns = ['A','C','G','T'], index = np.arange(1,c))
Sites = ['TAAT','ATTA','[W|A|T]{4}']


## Search for TAAT, ATTA, and flex AT rich sites
## Score motifs based on individual nucleotide selection of site

Found_sites_scores = []
nmatches = 0
for s in Sites:
    print('Starting search with', s)
    score = []
    Seq_Score = []
    for match in re.finditer(s,seq):
        nmatches +=1
        for m in np.arange(0,4):
            # re starts at 0 - need to account for that
            score.append(max(motif.loc[match.start()+1+m]))
        Seq_Score = round(sum(score)/len(score),3)
        Found_sites_scores.append([match.group(),match.start()+1, 
            match.end(), Seq_Score])
        
Found_sites_scores_df = pd.DataFrame(Found_sites_scores, 
    columns = ['Seq','Start','End','Score'], index = np.arange(1,nmatches+1))
Found_sites_scores_df = Found_sites_scores_df.drop_duplicates()
print(Found_sites_scores_df)

if nmatches < 2:
    print('Two sites not found')
    with open(export_path,'w') as log:
        log.write('N/A')
    exit()

## Find top site
top_site = Found_sites_scores_df.loc[Found_sites_scores_df.loc[:,'Score'].idxmax(axis = 'columns')]
Index = Found_sites_scores_df.loc[:,'Score'].idxmax(axis = 'columns')

# ## Remove top site from dataframe to find next max
Found_sites_scores_df = Found_sites_scores_df.drop(Index)
print(Found_sites_scores_df)

found = False
while found == False and Found_sites_scores_df.empty == False:
    top_site2 = Found_sites_scores_df.loc[Found_sites_scores_df.loc[:,'Score'].idxmax(axis = 'columns')]
    if ( top_site2['End'] in np.arange(top_site['Start']-1, top_site['End']+2) ) or ( top_site2['Start'] in np.arange(top_site['Start']-1, top_site['End']+2)):
        ## Drop overlap match and repeat
        Index = Found_sites_scores_df.loc[:,'Score'].idxmax(axis = 'columns')
        Found_sites_scores_df = Found_sites_scores_df.drop(Index)
        print(Found_sites_scores_df)
    else:
        found = True
    
print(top_site)
print(top_site2)

if top_site['Start'] < top_site2['Start']:
    ## top site is first site
    Spacer = top_site2['Start'] - top_site['End'] - 1
    dimer_site = top_site['Seq'] + str(Spacer) + 'N' + top_site2['Seq']
elif top_site['Start'] > top_site2['Start']:
    ## top2 site is first site
    Spacer = top_site['Start'] - top_site2['End'] - 1
    dimer_site = top_site2['Seq'] + str(Spacer) + 'N' + top_site['Seq']

if Spacer < 0:
    dimer_site = 'NA'

export_path = path + '/long_motif_consensus.txt'
with open(export_path,'w') as log:
    log.write(dimer_site)
