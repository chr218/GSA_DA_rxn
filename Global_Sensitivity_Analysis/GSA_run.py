# -*- coding: utf-8 -*-
"""
Created on Tue Jul 27 10:26:05 2021

@author: chr218

This program applies ANOVA
Global Sensitivity Analysis
for the [4+2] DA rxn between
ethene and isoprene.
"""
import numpy as np
import matplotlib.pyplot as plt
from gsMk import PCE
import sys
import numbers
import csv


def new_dat(fname):
    #Read data
    rows = []
    with open(fname,"r") as csvfile:
        csvreader = csv.reader(csvfile)
        for row in csvreader:
            rows.append(row[1:])
            

    tot_arr = np.asarray(rows,dtype=float)
    x = tot_arr[:,0:50]
    x = np.hstack((tot_arr[:,11:16],tot_arr[:,27:32],tot_arr[:,34:35],tot_arr[:,40:41],tot_arr[:,45:50]))
    y_DA = np.sum(np.hstack((tot_arr[:,44+18:48+18],tot_arr[:,56+18:59+18])),axis=1)#Sum catalytic steps where gas isoprene is consumed.
    
    return y_DA, x
    

y_DA, x = new_dat('output_0.5_shorter.csv')


#%%
#Train Data to PCE

r = 1.0 #range of perturbations. for example: 0.5 - (-0.5) = 1.0 [kJ/mol] 
ntarget = len(y_DA) #Total number of data
ntrain = int(np.floor(0.9*ntarget)) #Number of training points.-We chose 90% of data.
nvar = 17#Number of reaction steps/intermediates.
nord = 3#Highest order of PCE


y_DA = np.log10(y_DA)#log base 10 prevents overflow.

# shuffle the index
idx = np.arange(ntarget)
np.random.shuffle(idx)
x = x[idx, :]
y_DA = y_DA[idx]

# shuffle the index
#idx = np.arange(ntarget)
#np.random.shuffle(idx)
#x = x[idx, :]
#y_DA = y_DA[idx]

#Fit polynomial chaos expansion
pce = PCE(nvar=nvar, nord=nord, polytype='Legd', withdev=False)
pce.fit(x[:ntrain, :] / r, y_DA[:ntrain])#divide by "r" to keep [-1,1].

ypred = pce.predict(x[ntrain:, :] / r)

rmse = np.sqrt(np.mean((y_DA[ntrain:] - ypred) ** 2))

print('RMSE: ',rmse)#Print RMSE of PCE


#%% calculate sobol indices
ssum = 0
for i in range(nvar):
    st = pce.sobol_index([i])
    print('SU {0:2d}   {1:8.6f}'.format(i, st))
    ssum += st
print('Sum      {0:4.2f}'.format(ssum))

for i in range(nvar):
    for j in range(i+1, nvar):
        print('SU {0:2d}, {1:2d}   {2:8.6f}'.format(i, j, pce.sobol_index([i, j])))


sobol_ij = np.zeros([nvar, nvar])
for i in range(nvar):
    for j in range(nvar):
        if i == j:
            s = pce.sobol_index([i])
        else:
            s = pce.sobol_index([i, j])
        sobol_ij[i, j] = s
        sobol_ij[j, i] = s



# for i in range(20):
#     for j in range(i + 1, 20):
#         for k in range(j + 1, 20):
#             s = pce.sobol_index([i, j, k])
#             if s >= 1e-3:
#                 print('SU {0:2d}, {1:2d}, {2:2d}   {3:8.6f}'.format(i, j, k, s))




#%% Single Indices

sobol_i, soboltol_i = [], []
for i in range(nvar):
    s, st = pce.sobol_index([i]), pce.sobol_totalindex([i])
    sobol_i.append(s)
    soboltol_i.append(st)
    print('{0:5d}  {1:6.2f}  {2:6.2f}'.format(i, s, st))


