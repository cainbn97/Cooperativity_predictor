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
Fold_change_mon = []
Fold_change_dim = []


for c in np.arange(1,5):
    c = str(c)
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
                    Bg_percent_mon = float(Bg_percent_temp.split('%')[0])
                    
                    Fold_change = Target_percent_mon/Bg_percent_mon
                    Fold_change_mon.append(Fold_change)
                
                elif len(Consensus_seq) > 10:
                    ## Parse dimer
                    Consensus_seq_dim = Consensus_seq
                    Target_percent_temp = line.split('\t')[5]
                    Target_percent_dim = round(float(Target_percent_temp.split('%')[0]),2)
                
                    Bg_percent_temp = line.split('\t')[7]
                    Bg_percent_dim = float(Bg_percent_temp.split('%')[0])
                    
                    Fold_change = Target_percent_dim/Bg_percent_dim
                    Fold_change_dim.append(Fold_change)

## Fit fold change curves
fit_mon = np.polyfit(np.arange(1,5), Fold_change_mon, 1)
slope_mon = fit_mon[0]
yint_mon = fit_mon[1]

fit_dim = np.polyfit(np.arange(1,5), Fold_change_dim, 1)
slope_dim = fit_dim[0]
yint_dim = fit_dim[1]

## Assess cooperativity - NEEDS WORK
if slope_dim > slope_mon:
    print('Dimer had higher enrichment than monomer.')
    Cooperative = 'True'


## Plot curves
if Fold_change_mon[3] > Fold_change_dim[3]:
    Max = Fold_change_mon[3]
else:
    Max = Fold_change_dim[3]

cy = np.arange(1,5)
plt.plot(cy,Fold_change_mon,'ro', 
    cy, Fold_change_dim, 'bo',
    cy, cy*slope_mon + yint_mon, 'r:',
    cy, cy*slope_dim + yint_dim, 'b:')
M, = plt.plot(cy,Fold_change_mon,'ro')
D, = plt.plot(cy,Fold_change_dim,'bo')
plt.xlabel('Cycle')
plt.ylabel('Fold change')
plot_title = TF + ' SELEX Enrichment'
plt.title(plot_title)
plt.axis([0.5, 4.5, 0, Max + 1])
plt.legend([M, D] , ['Monomer', 'Dimer'])
filename = TF+'_Enrichment_plot.png'
plt.savefig(filename)

## Write to log
with open('/users/cainu5/SELEX_analysis/Run_summary.txt','a') as log:
    log.write(TF+'\t'+ Cooperative +
    '\t'+str(Consensus_seq_dim)+'\t'+str(round(slope_dim,2)) +
    '\t'+str(Consensus_seq_mon)+'\t'+str(round(slope_mon,2)) +
    '\t'+str(Fold_change_mon)+'\t'+str(Fold_change_dim)+'\n')

