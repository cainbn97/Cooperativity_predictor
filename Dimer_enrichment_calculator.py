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
    Fold_change_dim = [1]
    Bg_percent_mon = []; Bg_percent_dim = []
    Target_percent_mon = []
    Target_percent_dim = []
    notes = []
    
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
                    Consensus_seq = line.split('\t')[1]
                    nmatches = 0
                    for match in re.finditer('W|R|Y|S|K|M|B|D|H|V|N', Consensus_seq):
                        nmatches += 1
                    if nmatches > 1 :
                        Fold_change_mon.append(0)
                        Target_percent_mon.append(0)
                        continue
                    
                    ## Check lengths - some motifs will not be found at each cycle
                    if ( Cycle1 == True and len(Fold_change_mon) != int(c) ) or ( Cycle1 == False and len(Fold_change_mon)+1 != int(c) ): 
                        Fold_change_mon.append(0)
                        notes.append('Monomer motif missing on cycle '+str(int(c)-1))                                                           
                   
                    Target_percent_temp = line.split('\t')[5]
                    Target_percent_mon.append(round(float(Target_percent_temp.split('%')[0]),2))
                
                    Bg_percent_temp = line.split('\t')[7]
                    Bg_percent_mon.append(round(float(Bg_percent_temp.split('%')[0]),2))
                    
                    Fold_change = Target_percent_mon[-1]/Bg_percent_mon[0]
                    Fold_change_mon.append(Fold_change)
    return [Consensus_seq, Fold_change_mon, Target_percent_mon, Bg_percent_mon, notes]

def fit_lines(cy, Fold_change):
    def func(x, a):
        return a*x
    fit = curve_fit(func, cy, Fold_change)
    slope = fit[0]
    return slope[0]

def Rsq(cy, array, slope):
    RSS = 0; TSS = 0
    for i in range(0,len(array)):
            RSS += (array[i] - (cy[i]*slope))**2
            TSS += (array[i] - sum(array)/len(array))**2
    R2 = 1 - (RSS/TSS)
    return R2


## Prepare log file
Run_summary = path + '/' + TF + '_Enrichment_analysis_run_summary.txt'
with open(Run_summary,'w') as log:
    log.write('TF\tDimer site\tCooperativity factor\tConsensus dimer sequence'+
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
                log.write(TF+'\t'+'No long motif found'+'\n')
            print('\tNo long motif found. Exiting')
            exit()

d_count = 0

Cycle1 = True
dimer_motifs = path + '/dimer_[0-9].motif'
for dimers in sorted(glob.glob(dimer_motifs)):
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
                    Consensus_seq_temp = line.split('-')[1]
                    Consensus_seq = Consensus_seq_temp.split(',')[0]
                    
                    ## Check lengths - some motifs will not be found at each cycle
                    if ( Cycle1 == True and len(Fold_change_dim) != int(c) ) or ( Cycle1 == False and len(Fold_change_dim)+1 != int(c) ): 
                        Fold_change_dim.append(0)
                        notes.append('Dimer motif missing on cycle '+str(int(c)-1) )                                                       
                   
                    ## Parse dimer
                    Consensus_seq_dim = Consensus_seq
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
                    
    ## Check for possible oversaturation - remove fourth cycle               
    if ( Target_percent_mon[-2] > 50 or Target_percent_mon2[-2] > 50 ):
        notes.append('Cycle 4 half site fold change ('+str(round(Fold_change_mon[-1],2))+') masked to avoid saturation')
        notes.append('Cycle 4 half site fold change ('+str(round(Fold_change_mon2[-1],2))+') masked to avoid saturation')
        Fold_change_mon[-1] = 0
        Fold_change_mon2[-1] = 0

    if Target_percent_dim[-2] > 50:
        notes.append('Cycle 4 dimer fold change ('+str(round(Fold_change_dim[-1],2))+') masked to avoid saturation')
        Fold_change_dim[-1] = 0

    ## Read dimer sequence from consensus sequence search
    long_consensus_path = path + '/long_motif_consensus.txt'
    with open(long_consensus_path, 'r') as dimer:            
        dimer_site = dimer.readlines()[d_count]
        dimer_site = str(dimer_site).strip()

    ## Check for prevalence of dimer site at cycle 4 - where de novo motif found
    if Target_percent_dim[-1] < 5:
        Cooperative = 0
        notes.append('Dimer motif at Cycle 4 had an enrichment of <5%.')
        print('\tDimer prevalence less than 5% - exiting')
        export_path = 'top_dimer_kmer_motifs_' + dimers
        os.remove(export_path + '/motif1.jpwm')
        os.remove(export_path + '/motif2.jpwm')
        os.rmdir(export_path)
        os.remove(dimers + '.motif')
        with open(Run_summary,'a') as log:
            log.write(TF+'\t'+ str(dimer_site)+'\t'+str(Cooperative) + 
            '\t'+str(Consensus_seq_dim)+'\t\t\t'
            +str(Consensus_seq_mon)+'\t\t\t'+ str(np.round(Target_percent_mon))+'\t'
            +str(np.round(Target_percent_mon2)) + '\t' + 
            +str(np.round(Target_percent_dim))+'\t\t\t\t'+str(notes)+'\n')
            continue
                    
    ## Initialize arrays
    if Cycle1 == False:
        cy = np.arange(1,5)
    else:
        cy = np.arange(0,5)

    Fold_change_mon = np.array(Fold_change_mon)
    Fold_change_mon2 = np.array(Fold_change_mon2)
    Fold_change_dim = np.array(Fold_change_dim)
        
    
    ## Mask 0s (where homer did not find a result)
    mask_mon = (Fold_change_mon != 0)
    mask_mon2 = (Fold_change_mon2 != 0)
    mask_dim = (Fold_change_dim != 0)

    ## Natural log transform
    cy_mon = cy[mask_mon]
    Fold_change_mon = np.log(Fold_change_mon[mask_mon])
    
    cy_mon2 = cy[mask_mon2]
    Fold_change_mon2 = np.log(Fold_change_mon2[mask_mon2])

    cy_dim = cy[mask_dim]
    Fold_change_dim = np.log(Fold_change_dim[mask_dim])

    if np.all(Fold_change_mon == 0) and np.all(Fold_change_mon2[1:5] == 0 ):
        print('Monomer sites are too deprecated. Cooperativity factor calculation \
is not recommended. Exiting.')
        notes.append('Half site 1 ' + Consensus_seq_mon + ' is too ambiguous')
        notes.append('Half site 2 ' + Consensus_seq_mon2 + ' is too ambiguous')
        exit()
    elif np.all(Fold_change_mon == 0):
        print('\tHalf site 1 is too deprecated. Will not be used for slope calculation.')
        Fold_change_mon_avg = Fold_change_mon2
        notes.append('Half site 1  ' + Consensus_seq_mon + ' is too ambiguous')
    elif np.all(Fold_change_mon2[1:5] == 0 ):
        print('\tHalf site 2 is too deprecated. Will not be used for slope calculation.')
        Fold_change_mon_avg = Fold_change_mon
        notes.append('Half site 2 ' + Consensus_seq_mon2 + ' is too ambiguous')
    else:
        Fold_change_mon_avg = np.mean(np.array( [Fold_change_mon, Fold_change_mon2] \
            ), axis = 0)

    ## Fit lines   
    slope_mon = fit_lines(cy_mon, Fold_change_mon)
    slope_mon2 = fit_lines(cy_mon2, Fold_change_mon2)
    slope_mon_avg = fit_lines(cy_mon, Fold_change_mon_avg)
    slope_dim = fit_lines(cy_dim, Fold_change_dim)

    ## Calculate Rsquared
    R2_mon = Rsq(cy_mon, Fold_change_mon, slope_mon)
    R2_mon2 = Rsq(cy_mon2, Fold_change_mon2, slope_mon2)
    R2_dim = Rsq(cy_dim, Fold_change_dim, slope_dim)

    ## Assess cooperativity - NEEDS THRESHOLD
    if slope_dim > slope_mon_avg:
        print('\tDimer had higher enrichment than monomer.')

    Cooperative = np.round(slope_dim/slope_mon,2)

    ## Plot natural log transformed curves
    if Fold_change_mon[-1] > Fold_change_dim[-1]:
        Max = Fold_change_mon[-1]
    else:
        Max = Fold_change_dim[-1]
    
    plt.plot(cy_mon, Fold_change_mon,'ro', 
        cy_mon2, Fold_change_mon2, 'co',
        cy_dim, Fold_change_dim, 'bo',
        cy, cy*slope_mon_avg, 'm:',
        cy, cy*slope_dim, 'b:')
    M, = plt.plot(cy_mon, Fold_change_mon,'ro')
    M2,= plt.plot(cy_mon2, Fold_change_mon2, 'co')
    D, = plt.plot(cy_dim, Fold_change_dim,'bo')
    plt.xlabel('Cycle')
    plt.ylabel('ln(Fold change)')
    plot_title = TF + ' SELEX Enrichment of ' + dimer_site
    plt.title(plot_title)
    plt.axis([-0.25, 4.5, 0, Max + 1])
    
    if np.all(Fold_change_mon == 0):
        plt.legend([M, M2, D] , [f'Half site 1: N/A', \
            f'Half site 2: {np.round(slope_mon2,2)}', f'Dimer: {np.round(slope_dim,2)}'])
    elif np.all(Fold_change_mon2[1:5] == 0 ):
        plt.legend([M, M2, D] , [f'Half site 1: {np.round(slope_mon,2)}', \
            f'Half site 2: N/A', f'Dimer: {np.round(slope_dim,2)}'])
    else:
        plt.legend([M, M2, D] , [f'Half site 1: {np.round(slope_mon,2)}', \
            f'Half site 2: {np.round(slope_mon2,2)}', f'Dimer: {np.round(slope_dim,2)}'])
    
        
    filename = TF + '_NatLog_2_Enrichment_plot_' + dimers + '.png'
    plt.savefig(filename)
    plt.close()

    ## Write to log
    with open(Run_summary,'a') as log:
        log.write(TF+'\t'+ str(dimer_site)+'\t'+str(Cooperative) + 
        '\t'+str(Consensus_seq_dim)+'\t'+str(np.round(slope_dim,2))+'\t'+str(np.round(R2_dim,4))+
        '\t'+str(Consensus_seq_mon)+'\t'+str(np.round(slope_mon,2))+'\t'+str(np.round(R2_mon,4))+
        '\t'+str(Consensus_seq_mon2)+'\t'+str(np.round(slope_mon2,2))+'\t'+str(np.round(R2_mon2,4))+
        '\t'+str(np.round(Target_percent_mon,2))+'\t'+str(np.round(Target_percent_dim,2))+
        '\t'+str(np.round(Fold_change_mon,2))+'\t'+str(np.round(Fold_change_dim,2))+'\t'+str(notes)+'\n')
        
    d_count += 1