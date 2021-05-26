#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Dimer enrichment script
"""

import numpy as np
import os
import matplotlib.pyplot as plt
import scipy.optimize as opt
from scipy.interpolate import interp1d

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
Run_summary = '/users/cainu5/SELEX_analysis/Run_fitting_tests.txt'

## Check if homer run finished
homer_html = path+'/Cycle4/'+TF+'_'+'4_homer_denovo_long/homerResults.html'

if os.path.exists(homer_html) == False:
    with open(Run_summary,'a') as log:
        log.write(TF+'\tRun did not finish within 24 hours\n')
    print('Run did not complete.')
    exit()

# Check if long motif was found
file_top_long = path+'/Cycle4/'+TF+'_4_homer_denovo_long/top_long.txt'

if os.stat(file_top_long).st_size == 0:
    with open(Run_summary,'a') as log:
        log.write(TF+'\t'+'No long motif found'+'\n')
    print('No long motif found. Exiting')
    exit()

for c in np.arange(1,5):
    c = str(c)
    print('Starting cycle ', c)
    file = path+'/Cycle'+c+'/'+TF+'_'+c+'_homer/knownResults.txt'
       
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
                if len(Consensus_seq) < 9: 
                    ## Check lengths - some motifs will not be found at each cycle
                    if len(Fold_change_mon) != int(c):
                        Fold_change_mon.append(0)
                        notes.append('Monomer motif missing on cycle '+str(c))
                        
  
                    ## Parse monomer
                    Consensus_seq_mon = Consensus_seq
                    Target_percent_temp = line.split('\t')[6]
                    Target_percent_mon.append(round(float(Target_percent_temp.split('%')[0]),2))
                
                    Bg_percent_temp = line.split('\t')[8]
                    Bg_percent_mon.append(round(float(Bg_percent_temp.split('%')[0]),2))
                    
                    Fold_change = Target_percent_mon[-1]/Bg_percent_mon[0]
                    Fold_change_mon.append(Fold_change)
                    
                                        
                elif len(Consensus_seq) > 9:                           
                    ## Check lengths - some motifs will not be found at each cycle
                    if len(Fold_change_dim) != int(c):
                        Fold_change_dim.append(0)
                        notes.append('Monomer motif missing on cycle '+str(c))
                        
                    ## Parse dimer
                    Consensus_seq_dim = Consensus_seq
                    Target_percent_temp = line.split('\t')[6]
                    Target_percent_dim.append(round(float(Target_percent_temp.split('%')[0]),2))
                    print(Target_percent_dim)
                    Bg_percent_temp = line.split('\t')[8]
                    Bg_percent_dim.append(round(float(Bg_percent_temp.split('%')[0]),2))
                    
                    try:
                        Fold_change = Target_percent_dim[-1]/Bg_percent_dim[0]
                    except ZeroDivisionError:
                        Fold_change = Target_percent_dim[-1]/(Bg_percent_dim[0]+0.01)
                        notes.append('0.01 added to Cycle '+str(c) + 
                            ' to avoid divide by zero error.')
                    
                    Fold_change_dim.append(Fold_change)
                

## Check for possible oversaturation - remove fourth cycle               
# if Target_percent_mon[-2] > 50:
    # notes.append('Cycle 4 monomer fold change ('+str(round(Fold_change_mon[-1],2))+') masked to avoid saturation')
    # Fold_change_mon[-1] = 0

# if Target_percent_dim[-2] > 50:
    # notes.append('Cycle 4 dimer fold change ('+str(round(Fold_change_dim[-1],2))+') masked to avoid saturation')
    # Fold_change_dim[-1] = 0

## Check for prevalence of dimer site at cycle 3 - where de novo motif found
if Target_percent_dim[-2] < 1:
    ## Write to log
    dimer_site = 'N/A'
    CF = 0
    notes.append('Dimer motif at Cycle 3 had an enrichment of <1%.')
    print('Dimer prevalence less than 1% - exiting')
    with open(Run_summary,'a') as log:
        log.write(TF+'\t'+ str(dimer_site)+'\t'+str(CF) + 
        '\t'+str(Consensus_seq_dim)+'\t'+''+'\t'''+
        '\t'+str(Consensus_seq_mon)+'\t'+''+'\t'+''+'\t'+''+'\t'+str(notes)+
        '\n')
        exit()
                
## Initialize arrays
cy = np.arange(0,5)
Fold_change_mon = np.array(Fold_change_mon)
Fold_change_dim = np.array(Fold_change_dim)
print(Fold_change_mon); print(Fold_change_dim)

## Mask 0s (where homer did not find a result)
mask_mon = (Fold_change_mon != 0)
mask_dim = (Fold_change_dim != 0)

cy_mon = cy[mask_mon]
cy_dim = cy[mask_dim]

if max(Fold_change_mon) > max(Fold_change_dim):
    Max = max(Fold_change_mon)
else:
    Max = max(Fold_change_dim)

## Plot non-transformed data
plt.plot(cy_mon, Fold_change_mon[mask_mon],'ro', 
    cy_dim, Fold_change_dim[mask_dim], 'bo')
M, = plt.plot(cy_mon, Fold_change_mon,'ro')
D, = plt.plot(cy_dim, Fold_change_dim,'bo')
plt.xlabel('Cycle')
plt.ylabel('Fold change')
plot_title = TF + ' SELEX Enrichment'
plt.title(plot_title)
plt.axis([-0.25, 4.5, 0, Max + 50])
plt.legend([M, D] , ['Monomer', 'Dimer'])
filename = TF+'_Enrichment_plot.png'
plt.savefig(filename)


#######     Fit to LogIt Function       #########

## Interpolate more points so logit can be fit
cy_int = np.linspace(0,4,num = 100, endpoint = True)
mon_int = interp1d(cy_mon, Fold_change_mon, 'cubic')
dim_int = interp1d(cy_dim, Fold_change_dim, 'cubic')

## Define logit function
def f( x, L, k, x0):
    return L/(1 + np.exp(-k * (x-x0)))


(LM, kM, x0M), _ = opt.curve_fit(f, cy_int, mon_int(cy_int))
(LD, kD, x0D), _ = opt.curve_fit(f, cy_int, dim_int(cy_int), 
    p0 = [max(dim_int(cy_int)),1,1])

# Plot logIt function
fit_mon = f(cy_int, LM, kM, x0M)
fit_dim = f(cy_int, LD, kD, x0D)
M, = plt.plot(cy_mon, Fold_change_mon,'ro')
D, = plt.plot(cy_dim, Fold_change_dim,'bo')
M_fit, = plt.plot(cy_int, fit_mon, 'r:')
D_fit, = plt.plot(cy_int, fit_dim, 'b:')
plt.xlabel('Cycle')
plt.ylabel('Fold change')
plot_title = TF + ' SELEX Enrichment'
plt.title(plot_title)
plt.axis([-0.25, 4.5, 0, Max + 50])
plt.legend([M, D, M_fit, D_fit] , 
    ['Monomer', 'Dimer', f'Monomer Logistic Fit: k = {round(kM,3)}', f'Dimer Logistic Fit: k = {round(kD,3)}'])
filename = TF+'_LogIt_Enrichment_plot.png'
plt.savefig(filename)

CF = kD/kM

## Calculate Rsquared
RSS_mon = 0; TSS_mon = 0; RSS_dim = 0; TSS_dim = 0
for i in range(0,len(Fold_change_mon)):
    RSS_mon += (Fold_change_mon[i] - (f(cy_mon[i], LM, kM, x0M)))**2
    TSS_mon += (Fold_change_mon[i] - sum(Fold_change_mon)/len(Fold_change_mon))**2
for i in range(0,len(Fold_change_dim)):
    RSS_dim += (Fold_change_dim[i] - (f(cy_dim[i], LD, kD, x0D)))**2
    TSS_dim += (Fold_change_dim[i] - sum(Fold_change_dim)/len(Fold_change_dim))**2   
R2_mon = 1 - (RSS_mon/TSS_mon)
R2_dim = 1 - (RSS_dim/TSS_dim)

#####     Write to log          #####
# Read dimer sequence from consensus sequence search
long_consensus_path = path + '/long_motif_consensus.txt'
with open(long_consensus_path, 'r') as dimer:
    for line in dimer:
        dimer_site = line

## Write to log
with open(Run_summary,'a') as log:
    log.write(TF+'\t'+ str(dimer_site)+'\t'+ str(round(CF,3)) + 
    '\t'+str(Consensus_seq_dim)+'\t'+str(round(kD,2))+'\t'+str(round(R2_dim,4))+
    '\t'+str(Consensus_seq_mon)+'\t'+str(round(kM,2))+'\t'+str(round(R2_mon,4))+
    '\t'+str(np.round(Fold_change_mon,2))+'\t'+str(np.round(Fold_change_dim,2))+'\t'+str(notes)+'\n')