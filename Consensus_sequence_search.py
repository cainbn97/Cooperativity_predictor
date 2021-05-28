#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
cainu5
05/26/21

Detects top dimer site from Homer's de novo motif analysis
'''
import os
import numpy as np
import glob
import pandas as pd
import re
import argparse

## Determine run type
parser = argparse.ArgumentParser(description='Trim monomer motif')
parser.add_argument('-c', '--COSMO', action = 'store_true', default = False,
    required = False, help = 'Generate motifs for COSMO run')
args = parser.parse_args()
COSMO = bool(args.COSMO)

Site_length = 4
Top2SpacThres = 1.5

## Find and save necessary folders/paths
path = os.getcwd()
TF = os.path.basename(path)

## Determine if there is a dimer site present in any de novo motif analyses
D_site_found = False
de_novo_motif_folder = path + '/Cycle4/' + TF +'_4_homer_denovo_long/homerResults/motif[0-9].motif'
## Grab all motif files from de novo result
for de_novo_motifs in sorted(glob.glob(de_novo_motif_folder)):
    print('Starting ', os.path.basename(de_novo_motifs))
    with open(de_novo_motifs, 'r') as de_novo_motif_file:
        ## Read through motif file by line
        c = 0
        motif_final = []
        for line in de_novo_motif_file:
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
    motif = pd.DataFrame(motif_final, columns = ['A','C','G','T'], index = np.arange(1,c))

    ## Search top kmers
    ## Score motifs based on individual nucleotide selection of site

    Found_sites_scores = []
    Seq_Score = []
    
    ## Adjust frameshifts to grab overlapping motifs
    for k in np.arange(0,len(seq)+1- Site_length):
        score = []
        ## Search new frameshift for selected site
        for m in np.arange(k,k+Site_length):                    
            score.append(max(motif.loc[m+1,:]))
        
        ## Calculate average of nt strength - append to dataframe
        Seq_Score = round(sum(score)/len(score),3)
        Found_sites_scores.append([seq[k:k+Site_length],k+1, 
            Site_length+k, Seq_Score])

    Found_sites_scores_df = pd.DataFrame(Found_sites_scores, 
        columns = ['Seq','Start','End','Score'], index = np.arange(0,len(seq)+1-Site_length))
    Found_sites_scores_df = Found_sites_scores_df.drop_duplicates()

    ## Find top site
    top_site = Found_sites_scores_df.loc[Found_sites_scores_df.loc[:,'Score'].idxmax(axis = 'columns')]
    Index = Found_sites_scores_df.loc[:,'Score'].idxmax(axis = 'columns')

    ## Remove top site from dataframe to find next max
    Found_sites_scores_df = Found_sites_scores_df.drop(Index)
    Non_top_sites = Found_sites_scores_df

    found = False
    while found == False and Found_sites_scores_df.empty == False:
        top_site2 = Found_sites_scores_df.loc[Found_sites_scores_df.loc[:,'Score'].idxmax(axis = 'columns')]
        if ( top_site2['End'] in np.arange(top_site['Start']-1, top_site['End']+2) ) \
            or ( top_site2['Start'] in np.arange(top_site['Start']-1, top_site['End']+2)):
            ## Drop overlap match and repeat
            Index = Found_sites_scores_df.loc[:,'Score'].idxmax(axis = 'columns')
            Found_sites_scores_df = Found_sites_scores_df.drop(Index)
        else:
            found = True
        
    ## Compare scores from top kmers to scores of other bps
    Max_bp = motif.max(axis = 1)
    Max_bp = Max_bp.drop(range((top_site['Start']),(top_site['End']+1)))
    Max_bp = Max_bp.drop(range((top_site2['Start']),(top_site2['End']+1)))
    
    Top_sites_score = (top_site['Score'] + top_site2['Score'])/2
    Other_scores = Max_bp.sum(axis = 0) / Max_bp.shape[0]

    if Top_sites_score > Other_scores * Top2SpacThres:
        D_site_found = True
        export_path = path + '/Cycle4/GSX2_4_homer_denovo_long/D_site_motif.txt'
        with open(export_path, 'w') as log:
            log.write(os.path.basename(de_novo_motifs))
        break

## Define consensus sequence        
if top_site['Start'] < top_site2['Start']:
    ## top site is first site
    Spacer = top_site2['Start'] - top_site['End'] - 1
    dimer_site = top_site['Seq'] + str(Spacer) + 'N' + top_site2['Seq']
    motif1 = top_site
    motif2 = top_site2
elif top_site['Start'] > top_site2['Start']:
    ## top2 site is first site
    Spacer = top_site['Start'] - top_site2['End'] - 1
    dimer_site = top_site2['Seq'] + str(Spacer) + 'N' + top_site['Seq']
    motif1 = top_site2
    motif2 = top_site

export_path = path + '/long_motif_consensus.txt'
if D_site_found == False:
    print('No dimer sites found in any de novo motif analyses. Exiting')
    with open(export_path,'w') as log:
        log.write('N/A')
else:
    print(top_site); print(top_site2)
    with open(export_path,'w') as log:
        log.write(dimer_site)

## Generate motif files for COSMO        
if ( COSMO == True and D_site_found == True ):
    ## Index top sites from motif file
    motif1 = motif.loc[motif1['Start']:motif1['End']]*100
    motif1 = motif1.transpose()

    motif2 = motif.loc[motif2['Start']:motif2['End']]*100
    motif2 = motif2.transpose()

    os.mkdir('top_dimer_kmer_motifs')
    export_path_COSMO = path + '/top_dimer_kmer_motifs/'
    with open(export_path_COSMO + 'motif1.jpwm','w') as log:
        log.write(motif1.to_string(index = False, header = False))
    with open(export_path_COSMO + 'motif2.jpwm','w') as log:
        log.write(motif2.to_string(index = False, header = False))

