#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""
cainu5

Cooperativity quantification script

"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d

def plot_and_calc(i,c): 
    ## Gather all information
    print('\n\nEnter excel file with EMSA quantitation: ')
    Excel = input()
    print('\nEnter number of concentrations to parse: ')
    n_conc = int(input())
    print('\nEnter name of curve:  ')
    TF = str(input())
    print('\nEnter row number that table starts on: ')
    R = int(input())

    DF = pd.read_excel(Excel, sheet_name = 'Quantitation', usecols = 'A:D',
        skiprows = R-1, nrows = n_conc) 
    print(DF)

    Exp, = plt.plot(DF['Concentration'], DF['Average Number of Sites Filled'],c+'o',
        linestyle = 'None')

    ## Solve for Kd
    Conc_int = np.linspace(DF.loc[0,'Concentration'], DF.loc[ n_conc -1, 'Concentration'], num = 1000)
    SiteFilled = interp1d(DF['Concentration'], DF['Average Number of Sites Filled'])

    SiteFilled = pd.DataFrame({'Conc':np.round(Conc_int,2), 'AvgSites':np.round(SiteFilled(Conc_int),2)})
    Avg1 = SiteFilled.loc[abs(SiteFilled['AvgSites'] - 1).idxmin()]
    Kd = Avg1['Conc'].mean()
    print(Avg1); print(Kd)

    Conc = DF['Concentration'].to_numpy()
    AvgSite = DF['Average Number of Sites Filled'].to_numpy()

    ## Solve for C and f
    DistSq = 1; f = 1
    for C in np.arange(0,50.01,0.01):
        #for f in np.arange(0,1.01,0.01):
        a = Conc/Kd
        P1 = ((f*2*a)/(1+2*a+C*a**2))+((1-f)*(2*a))/(1+2*a)
        P2 = (f*C*a**2)/(1+2*a+C*a**2)
        Sitesf = (P1+P2*2)
        if sum((Sitesf - AvgSite)**2) < DistSq:
            DistSq = sum((Sitesf - AvgSite)**2)
            C_final = C; f_final = f; 
            Sitesf_final = Sitesf

    ## Print out outputs
    print(f'Cooperativity factor: {round(C_final,2)}')
    print(f'f factor: {round(f_final,2)}')
    print(f'Kd: {round(Kd,2)}')
    print(f'Average Euclidian distance: {round((np.mean((Sitesf_final - AvgSite)**2)),5)}')

    ## Final plot
    a = Conc_int/Kd
    P1 = ((f_final*2*a)/(1+2*a+C_final*a**2))+((1-f_final)*(2*a))/(1+2*a)
    P2 = (f_final*C_final*a**2)/(1+2*a+C_final*a**2)
    Sitesf = (P1+P2*2)
    
    i, = plt.plot(Conc_int, Sitesf, c+':', label = f'{TF} Experimental\nFitted curve: \nC: {round(C_final,2)}\nf: {round(f_final,2)}')

print('Enter the name of the experiment: ')
Experiment = str(input())

print('Enter number of curves you wish to plot: ')
curves = int(input())

for i in np.arange(curves):
    var = 'Fit_' + str(i)
    colors = ['b','r','g','c','m','y']
    plot_and_calc(var, colors[i])
   
fig_title = Experiment + '_cooperativity_calculation.png'
plt.xlabel(u"Concentration (\u03bcM)"); plt.ylabel('Average Number of Sites Filled')
# plt.errorbar(DF['Concentration'],DF['Average Number of Sites Filled'], xerr = None,
    # yerr = DF['STDEV'])
Title = Experiment + ' EMSA quantitation'
plt.title(Experiment)
plt.legend()
plt.savefig(fig_title)