#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Heatmap plotter for COSMO output
cainu5
06/10/21

"""

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import glob
import os

path = os.getcwd()
TF = os.path.basename(path)
dimer_paths = path + '/dimer_[0-9].motif'
orientation = ['motif1_motif1', 'motif1_motif2', 'motif2_motif2','motif2_motif1']
arrangement = ['FF','RF','FR']
cycle = []

if len(sorted(glob.glob(path + '/Cycle1/'+ TF + '_[0-9]_dimer_[0-9]_homer/knownResults.txt'))) > 0:
    end_index = 5
else:
    end_index = 4
 
long_consensus = []
with open(path + '/long_motif_consensus.txt','r') as long_consensus_file:
    for line in long_consensus_file:
        long_consensus.append(line)      
c = 0    

for dimers in sorted(glob.glob(dimer_paths)):
    dimer = os.path.basename(dimers)
    dimer = dimer.split('.')[0]
    print(f'\tStarting {long_consensus[c]} ...')
    for site_type in orientation:
        for site in arrangement: 

            files = path + '/Cycle[0-9]/*' + dimer + '_homer/*' + site_type + '_' + site + '.tab'
            matrix = pd.DataFrame(0, index=np.arange(0, end_index), columns = np.arange(1,11))
            cycle = []
            
            ## Create matrix of counts
            for file in sorted(glob.glob(files)):
                cycle_temp = file.split('/')[-1]
                cycle.append(cycle_temp.split('_')[0])
                with open(file, 'r') as readfile:
                    for line in readfile: 
                        spacer = int(line.split('|')[-2])
                        if spacer in np.arange(1,11):
                            count = int(line.split('|')[-1])
                            matrix.loc[len(cycle)-1,spacer] += count
                
            ## Convert NaNs to zeros in case there are missing values
            matrix = matrix.fillna(0)
            
            ## Check if matrix is all zeros
            if matrix.equals(pd.DataFrame(0, index=np.arange(end_index), columns = np.arange(1,11))):
                continue
            
            
            ## heatmap plotting
            plt.clf()
            filename = TF + '_' + dimer + '_' + site_type + '_'+ site + '_cosmo_output.png'
            HeatMap = sns.heatmap(matrix, annot = True, fmt = 'g')
            plt.xlabel('Spacer Length'); plt.ylabel('Cycle number')
            plt.title(TF + ' ' + str(long_consensus[c]) + ' '+ site + ' ' + ' Cosmo Output')
            fig = HeatMap.get_figure()
            # plt.figure(figsize=(12,6), dpi = 80)
            fig.savefig(filename)
            
            output = TF + '_' + dimer +  '_COSMO_counts_' + site_type + '_' + site + '.txt'
            with open(output, 'w') as log:
                log.write(matrix.to_string(index = True, header = True))
                
    c =+ 1