#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""

Dimer enrichment script
cainu5
06/10/21

"""

import numpy as np
import os
import matplotlib.pyplot as plt
import glob
from scipy.optimize import curve_fit
import re

## Grab working directory
path = os.getcwd()
TF = os.path.basename(path)

## Define functions
def parsing_half_sites(site):
    Cycle1 = True
    
    Fold_change_mon = [1]
    Bg_percent_mon = []
    Target_percent_mon = []
    Target_percent_dim = []
    notes = []
    
    
    ## Read consensus sequence directly from motif file    
    with open(path + '/PWMs_of_dimers/' + site + '_' + dimers + '.motif') as file:
        Consensus_seq = file.readline().split('\t')[1]
        nmatches = 0
        for match in re.finditer('W|R|Y|S|K|M|B|D|H|V|N', Consensus_seq):
            nmatches += 1
        if nmatches > 1 :
            Fold_change_mon = [1, 0, 0, 0, 0]
            Target_percent_mon = [0,0,0,0,0]
            print('\tHalf ' + site + ' was is too ambiguous...')
            notes.append('Half ' + site + ' was is too ambiguous...')
            return [Consensus_seq, Fold_change_mon, Target_percent_mon, Bg_percent_mon, notes]
    
    for c in np.arange(1,5):
        c = str(c)
        file = path+'/Cycle'+c+'/'+TF+'_'+c+'_' + str(site) + '_' + dimers + \
            '_mask_homer/knownResults.txt'
        ## Check if cycle 1 was used as background. Skip if it was
        if ( os.path.isfile(file) == False ) and ( c == '1'):
            Cycle1 = False
            continue
        
        ## Read in motif files
        with open(file) as readfile:
            for i, line in enumerate(readfile):
                ## Parse through folder
                if i == 0:
                    ## skip header
                    continue
                else:
                                                                                                
                    Target_percent_temp = line.split('\t')[5]
                    Target_percent_mon.append(round(float(Target_percent_temp.split('%')[0]),2))
                
                    Bg_percent_temp = line.split('\t')[7]
                    Bg_percent_mon.append(round(float(Bg_percent_temp.split('%')[0]),2))
                    
                    Fold_change = Target_percent_mon[-1]/Bg_percent_mon[0]
                    Fold_change_mon.append(Fold_change)
    
        ## Check lengths - some motifs will not be found at each cycle
        if ( Cycle1 == True and len(Fold_change_mon)-1 != int(c) ) or ( Cycle1 == False and len(Fold_change_mon) != int(c) ): 
            Fold_change_mon.append(0)
            Target_percent_mon.append(0)
            notes.append('Half ' + site + ' missing on cycle '+str(int(c))) 
            print('\tHalf ' + site + ' missing on cycle '+str(int(c)))
    
    if len(Fold_change_mon) < 3:
        print('\tHalf ' + site + ' was masked by dimer site...')
        notes.append('Half ' + site + ' was masked by dimer site.')
        Fold_change_mon = [1, 0, 0, 0, 0]
        Target_percent_mon = [0,0,0,0,0]
       
    return [Consensus_seq, Fold_change_mon, Target_percent_mon, Bg_percent_mon, notes]

def fit_lines(cy, Fold_change):
    def func(x, a):
        return a*x
    fit = curve_fit(func, cy, Fold_change)
    slope = fit[0]
    slope = slope[0]
    
    RSS = 0; TSS = 0
    for i in range(0,len(Fold_change)):
            RSS += (Fold_change[i] - (cy[i]*slope))**2
            TSS += (Fold_change[i] - sum(Fold_change)/len(Fold_change))**2
    R2 = 1 - (RSS/TSS)
    
    return slope, R2
    


## Prepare log file
Run_summary = path + '/' + TF + '_Enrichment_analysis_run_summary.txt'
with open(Run_summary,'w') as log:
    log.write('TF\tPWM\tDimer site\tCooperativity factor\tConsensus dimer sequence'+
    '\tDimer slope\tRsq\t' + 
    'Consensus half site 1\tHalf site 1 slope\tRsq\t'+
    'Consensus half site 2\tHalf site 2 slope\tRsq\t'+
    'Target percent half site 1\tTarget percent half site 2\t' + 
    'Fold change half site 1\tFold change half site 2\t'+
    'Fold change dimer\tNotes\n')

## Check if homer run finished
homer_html = path+'/Cycle4/'+TF+'_'+'4_homer_denovo_long/homerResults.html'

# Check if long motif was found
file_top_long = path + '/long_motif_consensus.txt'
with open(file_top_long,'r') as long_consensus:
    for line in long_consensus:
        if line == 'N/A':        
            with open(Run_summary,'a') as log:
                log.write(TF+'\t'+'No dimer site found'+'\n')
            print('\tNo dimer site found. Exiting')
            exit()

d_count = -1

Cycle1 = True
dimer_motifs = path + '/PWMs_of_dimers/dimer*'
for dimers in sorted(glob.glob(dimer_motifs)):
    d_count += 1
    
    print(f'\tStarting parsing for {os.path.basename(dimers)} homer runs...') 
    dimers = os.path.basename(dimers).split('.')[0]

    ## Initiate variables
    Fold_change_mon = [1]
    Fold_change_dim = [1]
    Bg_percent_mon = []; Bg_percent_dim = []
    Target_percent_mon = []
    Target_percent_dim = []
    notes = []
                        
    ## Parse half sites
    [Consensus_seq_mon, Fold_change_mon, Target_percent_mon, Bg_percent_mon1, notes_tmp] =\
        parsing_half_sites('site1')
    notes.append(notes_tmp)
    
    [Consensus_seq_mon2, Fold_change_mon2, Target_percent_mon2, Bg_percent_mon2, notes_tmp] =\
        parsing_half_sites('site2')
    notes.append(notes_tmp)

      
    for c in np.arange(1,5):
        c = str(c)
        file = path+'/Cycle'+c+'/'+TF+'_'+c+'_' + dimers + '_homer/knownResults.txt'
        
        with open(path + '/PWMs_of_dimers/' + dimers + '.motif') as dimer_motif:
            Consensus_seq_dim = dimer_motif.readline().split('\t')[1]
            Consensus_seq_dim = Consensus_seq_dim.split('-')[1]
            Consensus_seq_dim = Consensus_seq_dim.split(',')[0]
            
        ## Check if cycle 1 was used as background. Skip if it was
        if ( os.path.isfile(file) == False ) and ( c == '1'):
            Cycle1 = False
            continue
           
        ## Read in motif files
        with open(file) as readfile:
            for i, line in enumerate(readfile):
                ## Parse through folder
                if i == 0:
                    ## skip header
                    continue
                else:                                                    
                   
                    ## Parse dimer
                    Target_percent_temp = line.split('\t')[6]
                    Target_percent_dim.append(round(float(Target_percent_temp.split('%')[0]),2))
                    Bg_percent_temp = line.split('\t')[8]
                    Bg_percent_dim.append(round(float(Bg_percent_temp.split('%')[0]),2))
                    
                    try:
                        Fold_change = Target_percent_dim[-1]/Bg_percent_dim[0]
                    except ZeroDivisionError:
                        Fold_change = Target_percent_dim[-1]/(Bg_percent_dim[0]+0.01)
                        notes.append('0.01 added to Cycle '+str(int(c)-1) + 
                            ' to avoid divide by zero error.')
                    
                    Fold_change_dim.append(Fold_change)    
                    
        ## Check lengths - some motifs will not be found at each cycle
        if ( Cycle1 == True and len(Fold_change_dim)-1 != int(c) ) or ( Cycle1 == False and len(Fold_change_dim) != int(c) ): 
            Fold_change_dim.append(0)
            notes.append('Dimer motif missing on cycle '+str(int(c)) )
            print('\tDimer motif missing on cycle '+str(int(c)) )
    
    ## Read dimer sequence from consensus sequence search
    long_consensus_path = path + '/long_motif_consensus.txt'
    with open(long_consensus_path, 'r') as dimer:            
        dimer_site = dimer.readlines()[d_count]
        dimer_site = str(dimer_site).strip()

    ## Check for prevalence of dimer site at cycle 4 - where de novo motif found
    if len(Target_percent_dim) == 0:
        Cooperative = 0
        notes.append('Dimer motif at Cycle 4 had an enrichment of <5%.')
        print('\tDimer prevalence less than 5% - exiting')
        # export_path = 'top_dimer_kmer_motifs_' + dimers
        # os.remove(os.path.join(export_path, 'motif1.jpwm'))
        # os.remove(os.path.join(export_path, 'motif2.jpwm'))
        # os.rmdir(export_path)
        os.remove(path + '/PWMs_of_dimers/dimer_' + dimers + '.motif')
        with open(Run_summary,'a') as log:
            log.write(TF+'\t'+ str(os.path.basename(dimers)) + '\t'+ str(dimer_site)+'\t'+str(Cooperative) + 
            '\t'+str(Consensus_seq_dim)+'\t\t\t'
            +str(Consensus_seq_mon)+'\t\t\t'+ str(Consensus_seq_mon2) + 
            '\t\t\t' + str(np.round(Target_percent_mon,2))+'\t'
            +str(np.round(Target_percent_mon2,2)) + '\t'
            +str(np.round(Target_percent_dim,2))+'\t\t\t'+str(notes)+'\n')
        continue 
    
    
    if Target_percent_dim[-1] < 5:
        Cooperative = 0
        notes.append('Dimer motif at Cycle 4 had an enrichment of <5%.')
        print('\tDimer prevalence less than 5% - exiting')
        # export_path = 'top_dimer_kmer_motifs_' + dimers
        # os.remove(os.path.join(export_path, 'motif1.jpwm'))
        # os.remove(os.path.join(export_path, 'motif2.jpwm'))
        # os.rmdir(export_path)
        os.remove(path + '/PWMs_of_dimers/' + dimers + '.motif')
        with open(Run_summary,'a') as log:
            log.write(TF+'\t'+ str(os.path.basename(dimers)) + '\t'+ str(dimer_site)+'\t'+str(Cooperative) + 
            '\t'+str(Consensus_seq_dim)+'\t\t\t'
            +str(Consensus_seq_mon)+'\t\t\t'+ str(Consensus_seq_mon2) + 
            '\t\t\t' + str(np.round(Target_percent_mon,2))+'\t'
            +str(np.round(Target_percent_mon2,2)) + '\t'
            +str(np.round(Target_percent_dim,2))+'\t\t\t'+str(notes)+'\n')
        continue 
        
        
            
    
    if len(Target_percent_dim) < 3 or ( len(Target_percent_mon) < 3 and len(Target_percent_mon2) < 3):
        notes.append('Dimer appeared in ' + str(len(Target_percent_dim)) + ' cycles. Exiting')
        print('\tDimer appeared in ' + str(len(Target_percent_dim)) + ' cycles. Exiting')
        with open(Run_summary,'a') as log:
            log.write(TF+'\t'+ str(os.path.basename(dimers)) + '\t'+ str(dimer_site)+'\t' + str(0) + '\t\t\t\t\t\t\t\t\t\t\t\t\t\t' + 
                str(notes) + '\n')
        continue

    ## Initialize arrays
    if Cycle1 == False:
        cy = np.arange(1,5)
    else:
        cy = np.arange(0,5)

    ## Check for possible oversaturation - remove fourth cycle               
    if Target_percent_mon[-2] > 50:
        notes.append('Cycle 4 half site fold change ('+str(round(Fold_change_mon[-1],2))+') masked to avoid saturation')
        Fold_change_mon[-1] = 0
        
    if Target_percent_mon2[-2] > 50:
        notes.append('Cycle 4 half site fold change ('+str(round(Fold_change_mon2[-1],2))+') masked to avoid saturation')
        Fold_change_mon2[-1] = 0

    if Target_percent_dim[-2] > 50:
        notes.append('Cycle 4 dimer fold change ('+str(round(Fold_change_dim[-1],2))+') masked to avoid saturation')
        Fold_change_dim[-1] = 0
                    

    ## Natural log transform
    Fold_change_mon = np.log(Fold_change_mon).astype('float')
    Fold_change_mon[Fold_change_mon == np.NINF] = np.nan    
    
    Fold_change_mon2 = np.log(Fold_change_mon2).astype('float')
    Fold_change_mon2[Fold_change_mon2 == np.NINF] = np.nan  

    Fold_change_dim = np.log(Fold_change_dim).astype('float')
    Fold_change_dim[Fold_change_dim == np.NINF] = np.nan

    ## Neither half site is usable
    if np.all(np.isnan(Fold_change_mon[1:])) and np.all(np.isnan(Fold_change_mon2[1:])):
        print('\tExiting.')
        with open(Run_summary,'a') as log:
            log.write(TF+'\t'+ str(os.path.basename(dimers)) + '\t'+ str(dimer_site)+'\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t' + 
                str(notes) + '\n')
        continue
    
    Fold_change_mon_avg = np.nanmean(np.array( [Fold_change_mon, Fold_change_mon2] \
        ), axis = 0)
            
    ## Calculate slopes
    slope_mon, R2_mon = fit_lines(cy[~np.isnan(Fold_change_mon)], \
        Fold_change_mon[~np.isnan(Fold_change_mon)])
    slope_mon2, R2_mon2 = fit_lines(cy[~np.isnan(Fold_change_mon2)], \
        Fold_change_mon2[~np.isnan(Fold_change_mon2)])
    slope_mon_avg, R2_avg = fit_lines(cy[~np.isnan(Fold_change_mon_avg)], \
        Fold_change_mon_avg[~np.isnan(Fold_change_mon_avg)])        
    slope_dim, R2_dim = fit_lines(cy[~np.isnan(Fold_change_dim)], \
        Fold_change_dim[~np.isnan(Fold_change_dim)])
           
     
    
    ## Prep variables for logging
    slope_mon = str(np.round(slope_mon,2))
    slope_mon2 = str(np.round(slope_mon2,2))
    R2_mon = str(np.round(R2_mon,2))
    R2_mon2 = str(np.round(R2_mon2,2))
         
    ## Assess cooperativity - NEEDS THRESHOLD
    if slope_dim > slope_mon_avg:
        print('\tDimer had higher enrichment than monomer.')

    Cooperative = np.round(slope_dim/slope_mon_avg,2)

    ## Plot natural log transformed curves
    if Fold_change_mon_avg[~np.isnan(Fold_change_mon_avg)][-1] \
        > Fold_change_dim[~np.isnan(Fold_change_dim)][-1]:
        Max = Fold_change_mon_avg[~np.isnan(Fold_change_mon_avg)][-1]
    else:
        Max = Fold_change_dim[~np.isnan(Fold_change_dim)][-1]
    
    plt.plot(cy[~np.isnan(Fold_change_mon)], 
        Fold_change_mon[~np.isnan(Fold_change_mon)],'ro', 
        cy[~np.isnan(Fold_change_mon2)], 
        Fold_change_mon2[~np.isnan(Fold_change_mon2)], 'co',
        cy[~np.isnan(Fold_change_dim)], 
        Fold_change_dim[~np.isnan(Fold_change_dim)], 'bo',
        cy, cy*slope_mon_avg, 'm:',
        cy, cy*slope_dim, 'b:')
    plt.xlabel('Cycle')
    plt.ylabel('ln(Fold change)')
    plot_title = TF + ' SELEX Enrichment of ' + dimer_site
    plt.title(plot_title)
    plt.axis([-0.25, 4.5, 0, Max + 1])
    
    ## Set the legend
    M, = plt.plot(cy[~np.isnan(Fold_change_mon)], 
        Fold_change_mon[~np.isnan(Fold_change_mon)],'ro')
    M2,= plt.plot(cy[~np.isnan(Fold_change_mon2)], 
        Fold_change_mon2[~np.isnan(Fold_change_mon2)], 'co')
    D, = plt.plot(cy[~np.isnan(Fold_change_dim)], 
        Fold_change_dim[~np.isnan(Fold_change_dim)],'bo')
    
    if np.all(Fold_change_mon[1:5] == np.nan):
        plt.legend([M, M2, D] , [f'Half site 1: N/A', \
            f'Half site 2: {slope_mon2}', f'Dimer: {np.round(slope_dim,2)}'])
    elif np.all(Fold_change_mon2[1:5] == np.nan ):
        plt.legend([M, M2, D] , [f'Half site 1: {slope_mon}', \
            f'Half site 2: N/A', f'Dimer: {np.round(slope_dim,2)}'])
    else:
        plt.legend([M, M2, D] , [f'Half site 1: {slope_mon}', \
            f'Half site 2: {slope_mon2}', f'Dimer: {np.round(slope_dim,2)}'])
    
        
    filename = TF + '_NatLog_2_Enrichment_plot_' + dimers + '.png'
    plt.savefig(filename)
    plt.close()

    ## Write to log
    with open(Run_summary,'a') as log:
        log.write(TF+'\t'+ str(os.path.basename(dimers)) + '\t' + str(dimer_site)+'\t'+str(Cooperative) + 
        '\t'+str(Consensus_seq_dim)+'\t'+str(np.round(slope_dim,2))+'\t'+str(np.round(R2_dim,4))+
        '\t'+str(Consensus_seq_mon)+'\t'+slope_mon+'\t'+R2_mon+
        '\t'+str(Consensus_seq_mon2)+'\t'+slope_mon2+'\t'+R2_mon2+
        '\t'+str(np.round(Target_percent_mon,2))+'\t'+str(np.round(Target_percent_dim,2))+
        '\t'+str(np.round(Fold_change_mon,2))+'\t'+str(np.round(Fold_change_mon2,2))+'\t'+str(np.round(Fold_change_dim,2))+'\t'+str(notes)+'\n')
        
    print('\n\n')