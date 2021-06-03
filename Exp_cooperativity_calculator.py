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

print('Enter excel file with EMSA quantitation: ')
Excel = input()
print('Enter number of concentrations to parse: ')
n_conc = int(input())
TF = Excel.strip().split('_')[0]

DF = pd.read_excel(Excel, sheet_name = 'Cooperativity Factor', usecols = 'A:E',
    skiprows = 18, nrows = n_conc) 
print(DF)

Exp, = plt.plot(DF['Concentration'], DF['Average Number of Sites Filled'],'bo',
    linestyle = 'None')
plt.xlabel('Concentration (nM)'); plt.ylabel('Average Number of Sites Filled')
# plt.errorbar(DF['Concentration'],DF['Average Number of Sites Filled'], xerr = None,
    # yerr = DF['STDEV'])
Title = TF + ' EMSA quantitation'
plt.title(Title)

## Solve for Kd
Conc_int = np.linspace(DF.loc[0,'Concentration'], DF.loc[ n_conc -1, 'Concentration'], num = 1000)
SiteFilled = interp1d(DF['Concentration'], DF['Average Number of Sites Filled'])

SiteFilled = pd.DataFrame({'Conc':np.round(Conc_int,2), 'AvgSites':np.round(SiteFilled(Conc_int),2)})
Avg1 = SiteFilled.loc[abs(SiteFilled['AvgSites'] - 1).idxmin()]
Kd = Avg1['Conc'].mean()
print(Avg1); print(Kd)

Conc = DF['Concentration'].to_numpy()
AvgSite = DF['Average Number of Sites Filled'].to_numpy()

DistSq = 1
for C in np.arange(0,50.01,0.01):
    for f in np.arange(0,1.01,0.01):
        a = Conc/Kd
        P1 = ((f*2*a)/(1+2*a+C*a**2))+((1-f)*(2*a))/(1+2*a)
        P2 = (f*C*a**2)/(1+2*a+C*a**2)
        Sitesf = (P1+P2*2)
        if sum((Sitesf - AvgSite)**2) < DistSq:
            DistSq = sum((Sitesf - AvgSite)**2)
            C_final = C; f_final = f; 
            Sitesf_final = Sitesf;


print(f'Cooperativity factor: {round(C_final,2)}')
print(f'f factor: {round(f_final,2)}')
print(f'Kd: {round(Kd,2)}')
print(Sitesf_final); print(AvgSite)
print(f'Average Euclidian distance: {round((np.mean((Sitesf_final - AvgSite)**2)),5)}')

## Final plot
a = Conc_int/Kd
P1 = ((f_final*2*a)/(1+2*a+C_final*a**2))+((1-f_final)*(2*a))/(1+2*a)
P2 = (f_final*C_final*a**2)/(1+2*a+C_final*a**2)
Sitesf = (P1+P2*2)


fig_title = TF + '_cooperativity_calculation.png'
Fit, = plt.plot(Conc_int, Sitesf, 'b:')
plt.legend([Exp, Fit], ['Experimental', f'Fitted curve: \nC: {round(C_final,2)}\nf: {round(f_final,2)}'])
plt.savefig(fig_title)