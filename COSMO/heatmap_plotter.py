#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Heatmap plotter for COSMO output
cainu5
12/12/20
"""

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import glob
import os

path = os.getcwd()
TF = os.path.basename(path)
orientation = ['motif1_motif1', 'motif1_motif2', 'motif2_motif2','motif2_motif1']
arrangement = ['FF','RF','FR']
cycle = []

for site_type in orientation:
    for site in arrangement: 
        print('Starting ', site_type)
        files = '*/*' + site_type + '_' + site + '.tab'
        matrix = pd.DataFrame(0, index=np.arange(5), columns = np.arange(1,11))
        cycle = []
        
        ## Create matrix of counts
        for file in sorted(glob.glob(files)):
            print(file)
            cycle.append(file.split('/')[0])
            with open(file, 'r') as readfile:
                for line in readfile: 
                    spacer = int(line.split('|')[-2])
                    if spacer in np.arange(1,11):
                        count = int(line.split('|')[-1])
                        matrix.loc[len(cycle)-1,spacer] += count
            
        ## Convert NaNs to zeros in case there are missing values
        matrix = matrix.fillna(0)
        
        ## Check if matrix is all zeros
        if matrix.equals(pd.DataFrame(0, index=np.arange(5), columns = np.arange(1,11))):
            print('No count information at ', site_type,'.')
            continue
        
        
        ## heatmap plotting
        plt.clf()
        filename = TF + '_' + site_type + '_'+ site + '_cosmo_output.png'
        HeatMap = sns.heatmap(matrix, annot = True, fmt = 'g')
        plt.xlabel('Spacer Length'); plt.ylabel('Cycle number')
        plt.title(TF + ' ' + site_type + '_'+ site + ' Cosmo Output')
        fig = HeatMap.get_figure()
        fig.savefig(filename)
        
        output = path + '/' + TF + '_COSMO_counts_' + site_type + '_' + site + '.txt'
        with open(output, 'w') as log:
            log.write(matrix.to_string(index = True, header = True))