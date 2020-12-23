#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Chi-square goodness of fit test for COSMO output
What a disaster - this does not work whatsoever but scared to toss it

cainu5
12/21/20

"""

import os
from scipy.stats import chi2
import glob
import pandas as pd
import numpy as np
from scipy.stats import chisquare
from more_itertools import distinct_combinations
from scipy.stats import chi2_contingency

# def goodness_of_fit(array):
    # diff = []
    # chi = []
    # rownames = array.index
    # print(array.shape)
    
    # ## index through dataframe to grab differences between cycles
    # ## uses indices aka spacer length (1-10)
    # ## export an array - not indexed by spacer length
    # for i in np.arange(1,array.shape[1]+1):
        # diff.append(array.loc[rownames[1],i] - array.loc[rownames[0],i])
    
    # ## determine chi values - indices different between array and dataframe
    # for i in np.arange(0,array.shape[1]):
        # chi.append(round((diff[i]**2)/array.loc[rownames[0],i+1],2))
    # chi_tot = sum(chi)
    # p_value = chi2.sf(chi_tot,array.shape[1])
    # return(chi_tot, p_value)

path = os.getcwd()
TF = os.path.basename(path)
consensus = '/users/cainu5/SELEX_analysis/testing/' + TF + '/long_motif_consensus.txt'


## grab dimer consensus sequence
with open(consensus, 'r') as readfile:
    for line in readfile:
        dimer = line.strip()

site_type = 'motif1_motif1'
files = '*/*' + site_type + '_FF.tab'
matrix = pd.DataFrame(0, index=np.arange(5), columns = np.arange(1,11))
matrix_norm = pd.DataFrame(0, index=np.arange(5), columns = np.arange(1,11))
cycle = []

## Create matrix of counts
# for file in sorted(glob.glob(files)):
    # print(file)
    # #if 'Cycle0' in file or 'Cycle4' in file:
    # cycle = (file.split('/')[0])
    # print(cycle)
    # with open(file, 'r') as readfile:
        # for line in readfile: 
            # spacer = int(line.split('|')[-2])
            # if spacer in np.arange(1,11):
                # count = int(line.split('|')[-1])
                # matrix.loc[cycle,spacer] += count
                
                
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
print(matrix)

## Normalize both to proportions
rownames = matrix.index

for l in rownames:
    tot = sum(matrix.loc[l,:])
    print(tot)
    for i in np.arange(1,matrix.shape[1]+1):
        #matrix[l,:].sum
        matrix_norm.loc[l,i] = matrix.loc[l,i] / tot

print(matrix_norm)

## compare cycle 4 frequencies to cycle 0 frequencies
# chisq, p = chisquare(matrix_norm.loc['Cycle4'], matrix_norm.loc['Cycle0'],ddof = len(np.arange(1,11))-1)
# print(chisq); print(p)
# p_value = chi2.sf(chisq,matrix.shape[1])
# print(p_value)

chi2, p_ind, dof, ex = chi2_contingency(matrix, correction = True)
print(chi2, p_ind, dof)

all_combos = list(distinct_combinations(matrix_norm.columns,2))
print(all_combos)
for comb in all_combos:
    new_matrix = matrix_norm.loc[:,[comb[0],comb[1]]]
    #print(new_matrix)
    chi2, p, dof , ex = chi2_contingency(new_matrix)
    #p_value = chi2.sf(chisq,new_matrix.shape[1])
    print(f"Chi2 result for pair {comb}: {round(chi2,2)}, p-value: {round(p,2)}")
