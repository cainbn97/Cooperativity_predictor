#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Performs monomer site trimming for COSMO based on length and chosen setting
Note: if want motif shorter than 4bp select --top flag
'''
import os
import numpy as np
import glob
import pandas as pd
import re
import argparse
import math

## parse arguments from command line
parser = argparse.ArgumentParser(description='Trim monomer motif')
parser.add_argument('-ml','--mon_length', required = False, default = 4,
    help='Give desired motif length')
parser.add_argument('-l', '--left', action = 'store_true', default = False,
    required = False, help = 'Add flank to left of motif')
parser.add_argument('-r','--right', action = 'store_true', default = False,
    required = False, help = 'Add flank to right of motif')
parser.add_argument('-c', '--center', action = 'store_true', default = False,
    required = False, help = 'Centralize motif within flanks')
parser.add_argument('-t', '--top', action = 'store_true', default = False,
    required = False, help = 'Select top motif')

args = parser.parse_args()
mon_length = int(args.mon_length)
left = bool(args.left)
right = bool(args.right)
center = bool(args.center)
top = bool(args.top)

## Check arguments
if sum([left, right, center, top]) > 1:
    print('Please only enter one flag for flank location')
    exit()
elif sum([left, right, center, top]) == 0:
    top = True

if top == True:
    kmer_length = mon_length
else:
    kmer_length = 4

## Grab motif file
c = 0
motif_final = []
path = os.getcwd()
motif_file = path + '/*_short.motif'

## Read in motif
for file in glob.glob(motif_file):
    with open(file,'r') as readfile:
        for line in readfile:
            if c == 0:
                seq_temp = line.split('>')[1]
                seq = seq_temp.split('\t')[0]
                print(seq)
            else:
                motif_row_final = []
                motif_row = line.replace('\n','')
                for m in np.arange(0,4):
                    motif_row_final.append(float(motif_row.split('\t')[m]))
                motif_final.append(motif_row_final)          
            c = c + 1

if mon_length > len(seq):
    print('Please enter a length that is less than the monomer motif length (', 
        len(seq), 'bp)')
    exit()
         
motif = pd.DataFrame(motif_final, columns = ['A','C','G','T'], index = np.arange(1,c))

## Score motifs based on individual nucleotide selection of site

Found_sites_scores = []
nmatches = 0
seq_init = seq

Seq_Score = []
cut_nuc = 0

## Adjust frameshifts to grab overlapping motifs
while Found_sites_scores == [] and kmer_length > 2:
    cut_nuc = 0
    nmatches = 0
    Sites = '[W|A|T]{'+str(kmer_length)+'}'
    while len(seq) > kmer_length:
        seq = seq_init[cut_nuc:]
        ## Search new frameshift for selected site
        for match in re.finditer(Sites,seq):
            nmatches +=1
            score = []
            
            ## Score each AT rich site
            for m in np.arange(0,kmer_length):
                
                # restarts at 0 - need to account for that
                score.append(max(motif.loc[match.start()+m+cut_nuc+1,:]))
            
            ## Calculate average of nt strength - append to dataframe
            Seq_Score = round(sum(score)/len(score),3)
            Found_sites_scores.append([match.group(),match.start()+cut_nuc+1, 
                match.end()+cut_nuc, Seq_Score])
        cut_nuc += 1
    kmer_length -=1

if kmer_length < mon_length:
    print('Monomer sites shortened to ',kmer_length+1, 'find top AT rich sequence')

## Find top site
Found_sites_scores_df = pd.DataFrame(Found_sites_scores, 
    columns = ['Seq','Start','End','Score'], index = np.arange(1,nmatches+1))
Found_sites_scores_df = Found_sites_scores_df.drop_duplicates()
print(Found_sites_scores_df)
top_site = Found_sites_scores_df.loc[Found_sites_scores_df.loc[:,'Score'].idxmax(axis = 'columns')]  

## Calculate flanking lengths - if top flag not selected
if left == True:
    l_flank = mon_length - 4
    r_flank = 0
elif right == True: 
    l_flank = 0
    r_flank = mon_length - 4
elif center == True:
    l_flank = math.floor((mon_length - 4)/2)
    r_flank = math.ceil((mon_length - 4)/2)
else:
    ## Not important for top option
    l_flank = 0
    r_flank = 0

## Quick length check
if (int(top_site['Start']) - l_flank) < 0: 
    l_flank = int(top_site['Start'])
    r_flank = mon_length - 4 - l_flank
    print('Flank had to be changed due to length of motif')
    print('Left flank was changed to ', l_flank)
    print('Right flank was changed to ', r_flank)

if int(top_site['End']) + r_flank > len(seq_init):
    r_flank = len(seq_init) - int(top_site['End'])
    l_flank = mon_length - 4 - r_flank
    print('Flanks had to be changed due to length of motif')
    print('Left flank was changed to ',l_flank)
    print('Right flank was changed to ',r_flank)

motif1 = motif.loc[(top_site['Start']-l_flank):(top_site['End']+r_flank)]*100
motif1 = motif1.transpose()
print(motif1)

TF = os.path.basename(path)
export_path_COSMO = '/users/cainu5/SELEX_analysis/COSMO_output/' + TF + '/motifs/'
with open(export_path_COSMO + 'monomer_motif.jpwm','w') as log:
    log.write(motif1.to_string(index = False, header = False))
