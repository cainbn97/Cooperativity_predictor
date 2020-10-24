#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Dimer enrichment script
"""

import numpy as np
import os
import matplotlib.pyplot as plt
import math

## Grab working directory
path = os.getcwd()
TF = os.path.basename(path)

## Initiate variables
Fold_change_mon = []
Fold_change_dim = []

for c in np.arange(1,5):
    c = str(c)
    print('Starting cycle ', c)
    file = path+'/Cycle'+c+'/'+TF+'_'+c+'_homer/knownResults.txt'
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
                    ## Parse monomer
                    Consensus_seq_mon = Consensus_seq
                    Target_percent_temp = line.split('\t')[5]
                    Target_percent_mon = round(float(Target_percent_temp.split('%')[0]),2)
                
                    Bg_percent_temp = line.split('\t')[7]
                    Bg_percent_mon = round(float(Bg_percent_temp.split('%')[0]),2)
                    
                    Fold_change = Target_percent_mon/Bg_percent_mon
                    Fold_change_mon.append(Fold_change)
                    
                    ## Check lengths - some motifs will not be found at each cycle
                    if len(Fold_change_mon) != int(c):
                        Fold_change_mon.append(0)
                
                elif len(Consensus_seq) > 10:
                    ## Parse dimer
                    Consensus_seq_dim = Consensus_seq
                    Target_percent_temp = line.split('\t')[5]
                    Target_percent_dim = round(float(Target_percent_temp.split('%')[0]),2)
                
                    Bg_percent_temp = line.split('\t')[7]
                    Bg_percent_dim = round(float(Bg_percent_temp.split('%')[0]),2)
                    
                    Fold_change = Target_percent_dim/Bg_percent_dim
                    Fold_change_dim.append(Fold_change)
                    
                    ## Check lengths - some motifs will not be found at each cycle
                    if len(Fold_change_dim) != int(c):
                        Fold_change_dim.append(0)
                    

## Logarthmic transformation
cycles = np.arange(1,5)
log_fc_mon = np.log2(Fold_change_mon)
log_fc_dim = np.log2(Fold_change_dim)
cy = np.log2(cycles)
print(cycles)
print(log_fc_mon)
print(log_fc_dim)
print(cy)

## Fit log curves curves
fit_mon = np.polyfit(cy, log_fc_mon, 1)
slope_mon = fit_mon[0]
print(slope_mon)
yint_mon = fit_mon[1]
print(yint_mon)

fit_dim = np.polyfit(cy, log_fc_dim, 1)
slope_dim = fit_dim[0]
print(slope_dim)
yint_dim = fit_dim[1]
print(yint_dim)

## Assess cooperativity - NEEDS THRESHOLD
if slope_dim > slope_mon:
    print('Dimer had higher enrichment than monomer.')

Cooperative = round(slope_dim/slope_mon,2)

## Plot log transformed curves
if Fold_change_mon[3] > Fold_change_dim[3]:
    Max = log_fc_mon[3]
else:
    Max = log_fc_dim[3]

plt.plot(cy,log_fc_mon,'ro', 
    cy, log_fc_dim, 'bo',
    cy, cy*slope_mon+ yint_mon, 'r:',
    cy, cy*slope_dim + yint_dim, 'b:')
M, = plt.plot(cy,log_fc_mon,'ro')
D, = plt.plot(cy,log_fc_dim,'bo')
plt.xlabel('log2(Cycle)')
plt.ylabel('log2(Fold change)')
plot_title = TF + ' SELEX Enrichment'
plt.title(plot_title)
plt.axis([-0.25, 2.25, 0, Max + 1])
plt.legend([M, D] , ['Monomer', 'Dimer'])
filename = TF+'_Enrichment_plot.png'
plt.savefig(filename)

## Write to log
with open('/users/cainu5/SELEX_analysis/Run_summary.txt','a') as log:
    log.write(TF+'\t'+ str(Cooperative) +
    '\t'+str(Consensus_seq_dim)+'\t'+str(round(slope_dim,2)) +
    '\t'+str(Consensus_seq_mon)+'\t'+str(round(slope_mon,2)) +
    '\t'+str(np.round(Fold_change_mon,2))+'\t'+str(np.round(Fold_change_dim,2))+'\n')

