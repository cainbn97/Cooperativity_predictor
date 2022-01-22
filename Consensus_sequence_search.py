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
parser = argparse.ArgumentParser(description='Define consensus sequence')
parser.add_argument('-c', '--COSMO', action = 'store_true', default = False,
    required = False, help = 'Generate motifs for COSMO run')
parser.add_argument('-p', '--PWM', action = 'store_true', default = False,
    required = False, help = 'Use PWMs as input')
args = parser.parse_args()
COSMO = bool(args.COSMO)
PWM = bool(args.PWM)

Site_length = 4
Top2SpacThres = 1.6

## Find and save necessary folders/paths
path = os.getcwd()
TF = os.path.basename(path)

## Remove past files if there are any from past runs
if os.path.isfile(path + '/long_motif_consensus.txt') == True:
    os.remove(path + '/long_motif_consensus.txt')
    
if os.path.isfile(path + '/dimer_description_check.txt') == True:
    os.remove(path + '/dimer_description_check.txt')

if os.path.isfile(path + '/Cycle4/' + TF +  '_4_homer_denovo_long/D_site_motif.txt') == True:
    os.remove(path + '/Cycle4/' + TF + '_4_homer_denovo_long/D_site_motif.txt')

COSMO_output = path + '/top_dimer_kmer_motifs_dimer_[0-9]/'
for export_path_COSMO in sorted(glob.glob(COSMO_output)):
    if os.path.isdir(export_path_COSMO) == True:
        os.remove(export_path_COSMO + 'motif1.jpwm')
        os.remove(export_path_COSMO + 'motif2.jpwm')
        os.rmdir(export_path_COSMO)

## Prepare log file
with open(path + '/dimer_description_check.txt', 'a') as log :
    log.write('TF\tPWM\tDimer site\tTop site: Spacer/flank\tTop site average\tSpacer/flank score\tTop site 1 score\tTop site 2 score\tlogP\tTarget percentage\t Background percentage\n')

## Determine if there is a dimer site present in any de novo motif analyses
D_site_found = False
de_novo_motif_folder = path + '/Cycle4/' + TF +'_4_homer_denovo_long/PWMs/'
## Grab all motif files from de novo result
count=0


for de_novo_motifs in sorted(os.listdir(de_novo_motif_folder), key=str.casefold):
    print('Starting ', de_novo_motifs)
    motif_number = str(de_novo_motifs).split('.')[0]

    with open(os.path.join(de_novo_motif_folder,de_novo_motifs), 'r') as de_novo_motif_file:
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
    
    if len(seq) < (2*Site_length+2):
        print('\n',de_novo_motifs,' is shorter than required.\nCurrent PWM length is ', len(seq), '. \nThe required PWM length is ',(2*Site_length+2), '.')
        continue
    
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
    print(Found_sites_scores_df)

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
    Ratio = np.around(Top_sites_score/Other_scores, decimals = 1)
    
    if ( PWM == False ):
        with open(path + '/Cycle4/' + TF +'_4_homer_denovo_long/homerResults.txt','r') as read:
            for line in read:
                if line.split('\t')[0] == str(count+1):
                    if line.split('\t')[1] == '*':
                        logP = line.split('\t')[3]
                        Target = line.split('\t')[4]
                        Bg = line.split('\t')[5]
                    else:
                        logP = line.split('\t')[2]
                        Target = line.split('\t')[3]
                        Bg = line.split('\t')[4]
    else:
        logP = 'N/A'
        Target = 'N/A'
        Bg = 'N/A'
        
    if Ratio >= Top2SpacThres and top_site['Score'] >= 0.6 and top_site2['Score'] >= 0.6 :
        D_site_found = True
        export_path = path + '/Cycle4/' + TF + '_4_homer_denovo_long/D_site_motif.txt'
        with open(export_path, 'a') as log:
            log.write(de_novo_motifs)
            log.write('\n')
            
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
        print(top_site); print(top_site2)
        with open(export_path,'a') as log:
            log.write(dimer_site)
            log.write('\n')

        # Generate motif files for COSMO        
        if ( COSMO == True ):
            ## Index top sites from motif file
            motif1 = motif.loc[motif1['Start']:motif1['End']]*100
            motif1 = motif1.transpose()

            motif2 = motif.loc[motif2['Start']:motif2['End']]*100
            motif2 = motif2.transpose()

            export_path_COSMO = path + '/top_dimer_kmer_motifs_dimer_' + de_novo_motifs.split('.')[0] + '/'
            os.mkdir(export_path_COSMO)
            with open(export_path_COSMO + 'motif1.jpwm','w') as log:
                log.write(motif1.to_string(index = False, header = False))
            with open(export_path_COSMO + 'motif2.jpwm','w') as log:
                log.write(motif2.to_string(index = False, header = False))
    
    else:
        dimer_site = 'N/A'
    
    with open(path + '/dimer_description_check.txt', 'a') as log :
        log.write(TF+ '\t' + de_novo_motifs + '\t' + dimer_site + '\t' + str(round(Top_sites_score/Other_scores,3)) + '\t' + str(round(Top_sites_score,3)) + '\t'
        + str(round(Other_scores,3)) + '\t' + str(round(top_site['Score'],3)) + '\t' + str(top_site2['Score']) + '\t' + logP + '\t' + Target + '\t' + Bg + '\n')
    
    
    count+=1
    print('\n\n')

export_path = path + '/long_motif_consensus.txt'
if D_site_found == False:
    print('\tNo dimer sites found in any de novo motif analyses. Exiting')
    with open(export_path,'w') as log:
        log.write('N/A')