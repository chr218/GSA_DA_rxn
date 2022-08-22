# -*- coding: utf-8 -*-
"""
Created on Mon Dec 28 13:42:14 2020

Used to collect the Shomate parameters for the Diels-Alder MKM.

In reference to Updates 04/26/2021 onward
within "Diels Alder part 4)" notebook.

@author: chr218
"""

import pandas as pd
import numpy as np

from scipy import constants

import sys


'''
Constants
'''

kB = constants.physical_constants['Boltzmann constant'][0]
h = constants.physical_constants['Planck constant'][0]
R_J = constants.physical_constants['molar gas constant'][0]
R = 0.08205746#[L*atm/K/mol]
Na = constants.physical_constants['Avogadro constant'][0]

'''
Ensure site balance is preserved, no negative concentrations,
no negative pressures. Enruse pressure is ~ 1 [atm]
'''

def sanity_check(F, shift,T,v,micromol,abserr,relerr,timestep):
    
    #tolerance for total pressure.
    tol = 0.01
    
    #Ensure pressures are not negative and total pressure is ~ 1 [atm]
    P = 0
    for i,p in enumerate(F[shift:]):
        if p < 0:
            print('Error: negative pressure found in species %s!' % (i+shift))
            sys.exit()
        else:
            P += p*R*T/v*micromol
    if np.absolute(1-P) > tol:
        print('Error: Total pressure differs from 1 [atm], %s' % (P))
        
        
    #Ensure ads concentrations are not negative and total coverage = 1.    
    for i,o in enumerate(F[0:shift]):
        if o < 0:
            print('Error: negative ads concentration found in species %s!' % (i))
            sys.exit()
    O = np.sum(F[0:shift])
    if O > 1.1:
        print('Error: Site balance differs from 1, %s' % (O))

    #Find smallest non-zero value in F.
    m = min(i for i in F if i > 0)
    if m < abserr and m != 0:
        print('Warning: Smallest F value encountered (%s,%s) is smaller than the solvers absolute tolerance (%s) at timestep (%s) [hr]. ' %(i,m,abserr,timestep))
        print('If a component of F is approximately below atol, the number of correct digits is not guaranteed.')
        for i,f in enumerate(F):
            print(i,f)
        sys.exit()
    

'''
Remove NaN from pandas dataframe. Find better way to do this!
'''
def rm_nan(arr):
    indices = []
    
    if arr.ndim < 2:
        for i,row in enumerate(arr):
            if pd.isnull(row):
                indices.append([i])
    else:
        for i,row in enumerate(arr):
            if (pd.isnull(row)).any():
                indices.append([i])        
                
    return np.delete(arr,indices,axis=0)

'''
Calculate thermo values using Shomate Parameters
'''
def Shomate_S(param,T):
    t = T/1000.0
    S = param[0]*np.log(t) + param[1]*t + param[2]*t*t/2 + param[3]*t*t*t/3 - param[4]/(2*t*t) + param[6]
    return S


def Shomate_H(param,T):  
    t = T/1000.0
    H = param[0]*t + param[1]*t**2/2 + param[2]*t**3/3 + param[3]*t**4/4 - param[4]/t + param[5]
    return H

'''
Collect Shomate Parameters
'''
# =============================================================================
# Import Data
# ============================================================================= 
def get_Shomate(fname):
    data = pd.read_excel(fname)
    
    #Collect Gas- Phase Shomate Parameters.
    gas_min_keys = rm_nan(data.values[1:,0])
    params_gas_list = rm_nan(data.values[1:,1:9])
    params_gas_arr = np.asarray(params_gas_list)
    params_gas = params_gas_arr.astype(np.float)

    #Collect Ads- Phase Shomate Parameters, non TST
    local_min_keys = rm_nan(data.values[1:,10])
    params_ads_HO_list = rm_nan(data.values[1:,11:19])  
    params_ads_HO_arr = np.asarray(params_ads_HO_list)
    params_ads_HO = params_ads_HO_arr.astype(np.float)
    
    params_ads_2D_list = rm_nan(data.values[1:,21:29])  
    params_ads_2D_arr = np.asarray(params_ads_2D_list)
    params_ads_2D = params_ads_2D_arr.astype(np.float)  
    
    #Collect Ads- Phase TST Shomate Parameters.
    TST_keys = rm_nan(data.values[1:,30])
    params_TST_HO_list = rm_nan(data.values[1:,31:39])  
    params_TST_HO_arr = np.asarray(params_TST_HO_list)
    params_TST_HO = params_TST_HO_arr.astype(np.float)
    
    params_TST_2D_list = rm_nan(data.values[1:,41:49])  
    params_TST_2D_arr = np.asarray(params_TST_2D_list)
    params_TST_2D = params_TST_2D_arr.astype(np.float)
    
    return params_gas,params_ads_HO,params_ads_2D,params_TST_HO,params_TST_2D,local_min_keys,TST_keys,gas_min_keys


'''
Calculate ads/des approximations.
'''
def ads_des_approx(key,thermo_gas,thermo_HO,thermo_2D,T,h_):
    G_HO = 0
    G_2D = 0
    
    init, fin = key.split('-')
    reactant_indices = init.strip('()').split('+')
    product_indices = fin.strip('()').split('+')
    
    for idx in reactant_indices:
        if 'g' in idx:
            G_HO += (thermo_gas[1][int(idx.split('g')[0])]-T*thermo_2D[2][int(product_indices[0])]/1000+(h_))
            G_2D += (thermo_gas[1][int(idx.split('g')[0])]-T*thermo_2D[2][int(product_indices[0])]/1000+(h_))
        else:
            G_HO += thermo_HO[0][int(idx)]+T*thermo_HO[2][int(idx)]/1000
            G_2D += thermo_2D[0][int(idx)]+T*thermo_2D[2][int(idx)]/1000
        
    return G_HO, G_2D
               
'''
Parse Key
'''
def parse_key(key):
    
    #Seperate initial and final states
    init, fin = key.split('-')
    reactant_indices = init.strip('()').split('+')
    product_indices = fin.strip('()').split('+')
    
    #Check if reaction step is in gas-phase. (all states must have 'g')
    gas_phase_reaction = 0
    for state in (reactant_indices+product_indices):
        if 'g' in state:
            gas_phase_reaction += 1
            
    if gas_phase_reaction == len(reactant_indices+product_indices):
        gas_phase_reaction = True
    else:
        gas_phase_reaction = False
            
    return reactant_indices, product_indices, gas_phase_reaction
    
'''
Calculate thermo-data
'''  

def calc_thermo(params_gas,params_ads_HO,params_ads_2D,T,epsilon):
    
    #Calculate Gas- Phase values.
    H_gas = []
    S_gas = []
    G_gas = []
    for i,param in enumerate(params_gas):
        H_gas.append(Shomate_H(param,T))
        S_gas.append(Shomate_S(param,T))
        G_gas.append(H_gas[-1]-T*S_gas[-1]/1000.0)

    #Collect Ads- Phase values.
    H_ads_HO = []
    S_ads_HO = []
    G_ads_HO = []        
    #print('HO approximation data collected-OK')
    for i,param in enumerate(params_ads_HO):
        if 0 not in param:
            H_ads_HO.append(Shomate_H(param,T))
            S_ads_HO.append(Shomate_S(param,T))
            G_ads_HO.append(H_ads_HO[-1]-T*S_ads_HO[-1]/1000.0+(epsilon[i]))
        else:
            H_ads_HO.append(Shomate_H(param,T))
            S_ads_HO.append(Shomate_S(param,T))
            G_ads_HO.append(H_ads_HO[-1]-T*S_ads_HO[-1]/1000.0)
        
    #print('2D approximation data collected-OK')
    H_ads_2D = []
    S_ads_2D = []
    G_ads_2D = []        
    for i,param in enumerate(params_ads_2D):
        if 0 not in param:
            H_ads_2D.append(Shomate_H(param,T))
            S_ads_2D.append(Shomate_S(param,T))
            G_ads_2D.append(H_ads_2D[-1]-T*S_ads_2D[-1]/1000.0+(epsilon[i]))
        else:
            H_ads_2D.append(Shomate_H(param,T))
            S_ads_2D.append(Shomate_S(param,T))
            G_ads_2D.append(H_ads_2D[-1]-T*S_ads_2D[-1]/1000.0)           
            

    thermo_gas = [G_gas,H_gas,S_gas]
    thermo_HO = [G_ads_HO,H_ads_HO,S_ads_HO]                   
    thermo_2D = [G_ads_2D,H_ads_2D,S_ads_2D]  
    
    return thermo_gas, thermo_HO, thermo_2D
 
'''
Calculate TST Equilibrium Constants and rate constants.
'''

def calc_step_constants(thermo, keys,T,micromol,hr):
    
    k_f = []
    k_r = []
    
    G_f = []
    G_r = []
    
    H_f = []
    H_r = []
    
    S_f = []
    S_r = []

    A = (kB*T/h)*hr #[hr^-1], Otherwise [s^-1]
    nu = (kB/R/h)/(micromol)*hr #[micromol/L/atm/hr] for gas-phase reactions. Otherwise [mol/L/atm/s]
    
    #Calculate rate constants
    for i, key in enumerate(keys[1]):
        reactant_indices, product_indices, gas_phase_reaction = parse_key(key)

        #Calculate foreward barrier del_G
        init_G = 0
        init_H = 0
        init_S = 0
        for idx in reactant_indices:
            if 'g' in idx:
                init_G += (thermo[0][0][int(idx.strip('g')[0])])
                init_H += (thermo[0][1][int(idx.strip('g')[0])])
                init_S += (thermo[0][2][int(idx.strip('g')[0])])
            else:
                init_G += (thermo[1][0][int(idx)])
                init_H += (thermo[1][1][int(idx)])
                init_S += (thermo[1][2][int(idx)])
                            
        #Gas-Phase reactions use'nu' instead of kBT/h!        
        if gas_phase_reaction:
            G_f.append((thermo[2][0][i]-init_G))
            H_f.append((thermo[2][1][i]-init_H))
            S_f.append((thermo[2][2][i]-init_S))
            k_f.append(nu*np.exp((thermo[2][0][i]-init_G)*(-1/(R_J/1000)/T)))
            #print(key,k_f[-1])
            
        else:
            G_f.append((thermo[2][0][i]-init_G))
            H_f.append((thermo[2][1][i]-init_H))
            S_f.append((thermo[2][2][i]-init_S))
            k_f.append(A*np.exp((thermo[2][0][i]-init_G)*(-1/(R_J/1000)/T)))
            #print(key,k_f[-1])
            

        #Calculate reverse barrier del_G
        fin_G = 0
        fin_H = 0
        fin_S = 0
        for idx in product_indices:
            if 'g' in idx:
                fin_G += (thermo[0][0][int(idx.strip('g')[0])])
                fin_H += (thermo[0][1][int(idx.strip('g')[0])])
                fin_S += (thermo[0][2][int(idx.strip('g')[0])])
            else:
                fin_G += (thermo[1][0][int(idx)])
                fin_H += (thermo[1][1][int(idx)])
                fin_S += (thermo[1][2][int(idx)])
                
        #Gas-Phase reactions use'nu' instead of kBT/h!
        if gas_phase_reaction:
            G_r.append((thermo[2][0][i]-fin_G))
            H_r.append((thermo[2][1][i]-fin_H))
            S_r.append((thermo[2][2][i]-fin_S))
            k_r.append(nu*np.exp((thermo[2][0][i]-fin_G)*(-1/(R_J/1000)/T)))
        else:
            G_r.append((thermo[2][0][i]-fin_G))
            H_r.append((thermo[2][1][i]-fin_H))
            S_r.append((thermo[2][2][i]-fin_S))
            k_r.append(A*np.exp((thermo[2][0][i]-fin_G)*(-1/(R_J/1000)/T)))

    return k_f, k_r, [G_f,G_r], [H_f,H_r], [S_f,S_r]
