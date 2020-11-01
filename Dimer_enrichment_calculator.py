#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Dimer enrichment script
"""

import numpy as np
import os
import matplotlib.pyplot as plt

## Grab working directory
path = os.getcwd()
TF = os.path.basename(path)

## Initiate variables
Fold_change_mon = [1]
Fold_change_dim = [1]
Bg_percent_mon = []; Bg_percent_dim = []
Target_percent_mon = []
Target_percent_dim = []
notes = []

for c in np.arange(1,5):
    c = str(c)
    print('Starting cycle ', c)
    file = path+'/Cycle'+c+'/'+TF+'_'+c+'_homer/knownResults.txt'
    
    ## Check if long motif was found
    if os.path.exists(file) == False:
        with open('/users/cainu5/SELEX_analysis/Run_summary_103120.txt','a') as log:
            log.write(TF+'\t'+'No long motif found'+'\n')
        print('No long motif found. Exiting')
        exit()
       
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
                if len(Consensus_seq) < 10:                   
                    ## Check lengths - some motifs will not be found at each cycle
                    if len(Fold_change_mon) != int(c):
                        Fold_change_mon.append(0)
                        notes.append('Monomer motif missing on cycle'+str(c))
                        
  
                    ## Parse monomer
                    Consensus_seq_mon = Consensus_seq
                    Target_percent_temp = line.split('\t')[6]
                    Target_percent_mon.append(round(float(Target_percent_temp.split('%')[0]),2))
                
                    Bg_percent_temp = line.split('\t')[8]
                    Bg_percent_mon.append(round(float(Bg_percent_temp.split('%')[0]),2))
                    
                    Fold_change = Target_percent_mon[-1]/Bg_percent_mon[0]
                    Fold_change_mon.append(Fold_change)
                    
                                        
                elif len(Consensus_seq) > 10:                    
                                        
                    ## Check lengths - some motifs will not be found at each cycle
                    if len(Fold_change_dim) != int(c):
                        Fold_change_dim.append(0)
                        notes.append('Monomer motif missing on cycle'+str(c))
                        
                    ## Parse dimer
                    Consensus_seq_dim = Consensus_seq
                    Target_percent_temp = line.split('\t')[6]
                    Target_percent_dim.append(round(float(Target_percent_temp.split('%')[0]),2))
                
                    Bg_percent_temp = line.split('\t')[8]
                    Bg_percent_dim.append(round(float(Bg_percent_temp.split('%')[0]),2))
                    
                    Fold_change = Target_percent_dim[-1]/Bg_percent_dim[0]
                    Fold_change_dim.append(Fold_change)
                

## Check for possible oversaturation - remove fourth cycle               
if Target_percent_mon[-2] > 50: 
    Target_percent_mon[-1] = 0

if Target_percent_dim[-2] > 50:
    Target_percent_dim[-1] = 0


## Check for prevalence of dimer site at cycle 3 - where de novo motif found
if Target_percent_dim[-2] < 1:
    ## Write to log
    dimer_site = 'N/A'
    Cooperative = 0
    notes.append('Dimer motif at Cycle 3 had an enrichment of <1%.')
    print('Dimer prevalence less than 1% - exiting')
    with open('/users/cainu5/SELEX_analysis/Run_summary_103120.txt','a') as log:
        log.write(TF+'\t'+ str(dimer_site)+'\t'+str(Cooperative) + 
        '\t'+str(Consensus_seq_dim)+'\t'+''+'\t'''+'\t'+''+
        '\t'+str(Consensus_seq_mon)+'\t'+''+'\t'+''+'\t'+''+'\t'+str(notes)+
        '\n')
        exit()
                

## Initialize arrays
cy = np.arange(0,5)
Fold_change_mon = np.array(Fold_change_mon)
Fold_change_dim = np.array(Fold_change_dim)

## Mask 0s (where homer did not find a result)
mask_mon = (Fold_change_mon != 0)
mask_dim = (Fold_change_dim != 0)

## Natural log transform
cy_mon = cy[mask_mon]
Fold_change_mon = np.log(Fold_change_mon[mask_mon])

cy_dim = cy[mask_dim]
Fold_change_dim = np.log(Fold_change_dim[mask_dim])

## Fit lines
fit_mon = np.polyfit(cy_mon, Fold_change_mon, 1)
slope_mon = fit_mon[0]
yint_mon = fit_mon[1]

fit_dim = np.polyfit(cy_dim, Fold_change_dim, 1)
slope_dim = fit_dim[0]
yint_dim = fit_dim[1]

## Calculate Rsquared
RSS_mon = 0; TSS_mon = 0; RSS_dim = 0; TSS_dim = 0
for i in range(0,len(Fold_change_mon)):
    RSS_mon += (Fold_change_mon[i] - (cy_mon[i]*slope_mon + yint_mon))**2
    TSS_mon += (Fold_change_mon[i] - sum(Fold_change_mon)/len(Fold_change_mon))**2
for i in range(0,len(Fold_change_dim)):
    RSS_dim += (Fold_change_dim[i] - (cy_dim[i]*slope_dim + yint_dim))**2
    TSS_dim += (Fold_change_dim[i] - sum(Fold_change_dim)/len(Fold_change_dim))**2   
R2_mon = 1 - (RSS_mon/TSS_mon)
R2_dim = 1 - (RSS_dim/TSS_dim)

## Assess cooperativity - NEEDS THRESHOLD
if slope_dim > slope_mon:
    print('Dimer had higher enrichment than monomer.')

Cooperative = round(slope_dim/slope_mon,2)

## Plot natural log transformed curves
if Fold_change_mon[-1] > Fold_change_dim[-1]:
    Max = Fold_change_mon[-1]
else:
    Max = Fold_change_dim[-1]

plt.plot(cy_mon, Fold_change_mon,'ro', 
    cy_dim, Fold_change_dim, 'bo',
    cy, cy*slope_mon+ yint_mon, 'r:',
    cy, cy*slope_dim + yint_dim, 'b:')
M, = plt.plot(cy_mon, Fold_change_mon,'ro')
D, = plt.plot(cy_dim, Fold_change_dim,'bo')
plt.xlabel('Cycle')
plt.ylabel('ln(Fold change)')
plot_title = TF + ' SELEX Enrichment'
plt.title(plot_title)
plt.axis([-0.25, 4.5, 0, Max + 1])
plt.legend([M, D] , ['Monomer', 'Dimer'])
filename = TF+'_NatLog_2_Enrichment_plot.png'
plt.savefig(filename)

## Read dimer sequence from consensus sequence search
long_consensus_path = path + '/long_motif_consensus.txt'
with open(long_consensus_path, 'r') as dimer:
    for line in dimer:
        dimer_site = line

## Write to log
with open('/users/cainu5/SELEX_analysis/Run_summary_103120.txt','a') as log:
    log.write(TF+'\t'+ str(dimer_site)+'\t'+str(Cooperative) + 
    '\t'+str(Consensus_seq_dim)+'\t'+str(round(slope_dim,2))+'\t'+str(round(R2_dim,4))+
    '\t'+str(Consensus_seq_mon)+'\t'+str(round(slope_mon,2))+'\t'+str(round(R2_mon,4))+
    '\t'+str(np.round(Fold_change_mon,2))+'\t'+str(np.round(Fold_change_dim,2))+'\t'+str(notes)+'\n')