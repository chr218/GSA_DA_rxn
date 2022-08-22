# -*- coding: utf-8 -*-
"""
Created on Wed Nov 25 17:36:21 2020
NOTE: all imported values are referenced w.r.t gas-phase reactants.
@author: chr218
"""
from data_config import get_Shomate, calc_thermo, calc_step_constants, ads_des_approx, parse_key, sanity_check
from scipy.integrate import solve_ivp
from scipy import constants

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import time
import sys

t0 = time.time()
'''
Constants
'''

kB = constants.physical_constants['Boltzmann constant'][0]
h = constants.physical_constants['Planck constant'][0]
R_J = constants.physical_constants['molar gas constant'][0]
R = 0.08205746#[L*atm/K/mol]
Na = constants.physical_constants['Avogadro constant'][0]

'''
Functions & Subroutines
'''

def CSTR_MKM(h_,h_intermediates):
    
    #Collect rates-find better way to do this!
    rates = []
    rates_f = []
    rates_r = []
    def vectorfield(t,F,k,keys):
        '''
        Arguments:
            F:  vector of the gas & ads phase species.
            t: time
            kf:  vector of the ads-phase rate parameters (foreword)
            kr:  vector of the ads-phase rate parameters (reverse)
            
            F[0:24] are ads species, F0[25:34] are gas species.
            
        '''
    
        #Initialize rate list.-Corresponds to TST_keys.
        r_f = []
        r_r = []
    
        #Initialize mass-balance list.
        dFdt = []       
        
        #Define rate equation for each step.
        for i,key in enumerate(TST_keys):
            reactant_indices, product_indices, gas_phase_reaction = parse_key(key)
    
            #Define rates
            #Foreward
            r_f.append(k[0][i])
            for idx in reactant_indices:
                if 'g' in idx:
                    r_f[-1] *= F[int(idx.split('g')[0])+shift]*R*T/v*micromol #Make pressure (R*T/v)
                else:
                    r_f[-1] *= F[int(idx)]
    
            #Reverse
            r_r.append(k[1][i])
            for idx in product_indices:
                if 'g' in idx:
                    r_r[-1] *= F[int(idx.split('g')[0])+shift]*R*T/v*micromol #Make pressure (RT/v)
                else:
                    r_r[-1] *= F[int(idx)]
             
        #Rate of each step is (foreward - reverse). We convert lists to arrays and subtract.
        r_f = np.asarray(r_f)
        r_r = np.asarray(r_r)
        rates_f.append(r_f)
        rates_r.append(r_r)
        
        r_arr = r_f - r_r
        rates.append(r_arr)
        
        #Define Mass Balance ODEs
        #Ads Phase Species
        for i,species in enumerate(local_min_keys):
    
            r_ads = 0
            for j,step in enumerate(TST_keys):
                reactant_indices, product_indices, gas_phase_reaction = parse_key(step)
                        
                #Skip gas-phase reactions        
                if gas_phase_reaction:
                    continue
    
                for idx in reactant_indices:
                    if str(species) == idx:
                        r_ads -= r_arr[j]
                        
                for idx in product_indices:
                    if str(species) == idx:
                        r_ads += r_arr[j]
                        
            dFdt.append(r_ads)#No inlet or effluent because surface species.
    
        #Gas Phase Species
        for i,species in enumerate(gas_min_keys):
    
            r_ads = 0
            r_gas = 0
            for j,step in enumerate(TST_keys):
                
                reactant_indices, product_indices, gas_phase_reaction = parse_key(step)
    
                if gas_phase_reaction:
                    if str(species) in reactant_indices:
                        r_gas -= r_arr[j]
                    elif str(species) in product_indices:
                        r_gas += r_arr[j] 
    
                else:
                    if str(species) in reactant_indices:
                        r_ads -= r_arr[j]
                    elif str(species) in product_indices:
                        r_ads += r_arr[j]                
    
            dFdt.append((F_in[i+shift]-F[i+shift]+W*Sigma*r_ads+V_g*r_gas)*(1/tau))#Only Gas-Phase species enter/leave reactor.
    
    
        dFdt.append((F_in[-1]-F[-1])*(1/tau))#Append inert.
           
        return dFdt
    
    '''
    Program Begins Here
    '''
    # =============================================================================
    # Control Variables
    # =============================================================================
    
    #Scaling Variables. Change [mol] to [micromol] and [s] to [hr]
    micromol = 1.0E-6
    hr = 60*60
    
    #micromol = 1
    #hr = 1
    
    T = 368.15#Temperature [K]
    #P_tot = 1.0 #Total Pressure [atm]
    v = 0.0005*hr #Flowrate [L/hr], otherwise [L/s]
    W = 0.1 #mass of catalyst [g]
    #Sigma = 0.000173371220724102/(micromol)#[micromol-acid-sites/g-cat], otherwise [mol-acid-sites/g-cat], 1 acid site/unit-cell.
    Sigma = 0.000346748332834016/(micromol)#[micromol-acid-sites/g-cat], otherwise [mol-acid-sites/g-cat], 47 Si/Al.
    
    #rho = 0.000559719876242191 # Volume of unit cell per mass in [L/g-cat], 1 acid-site/unit-cell.
    rho = 0.000559729358492405# Volume of unit cell per mass in [L/g-cat], 47 Si/Al.
    
    V_cat = W*rho #Volume of catalyst inreactor [L]
    V_g = V_cat #Assuming equivalent volume for gas-phase. [L]
    tau = V_g/v #Space velocity #[hr], otherwise [s]
    
    # =============================================================================
    # Import Data
    # =============================================================================
    
    #Collect Shomate Parameters.
    params_gas,params_ads_HO,params_ads_2D,params_TST_HO,params_TST_2D,local_min_keys,TST_keys,gas_min_keys = get_Shomate("Shomate_ordered.xlsx")
    
    #Calculate Thermo values from Shomate parameters.
    #NOTEL thermo_HO[0,1,2] = G,H,S
    thermo_gas, thermo_HO, thermo_2D = calc_thermo(params_gas, params_ads_HO, params_ads_2D,T,h_intermediates)
    thermo_gas, thermo_HO_TST, thermo_2D_TST = calc_thermo(params_gas, params_TST_HO, params_TST_2D,T,h_)
    
    #Calculate ads/des TST thermo using approximation.
    for i,x in enumerate(thermo_HO_TST[0]):
        if not x:
            thermo_HO_TST[0][i], thermo_2D_TST[0][i] = ads_des_approx(TST_keys[i], thermo_gas, thermo_HO, thermo_2D, T, h_[i])
    
    thermo_total_HO = [thermo_gas,thermo_HO,thermo_HO_TST]
    thermo_total_2D = [thermo_gas,thermo_2D,thermo_2D_TST]
    keys = [gas_min_keys, TST_keys, local_min_keys]
    
    #Calculate rate constant for each step.
    k_f_HO, k_r_HO, G_HO_list, H_HO_list ,S_HO_list   = calc_step_constants(thermo_total_HO,keys,T,micromol,hr)
    k_f_2D, k_r_2D, G_2D_list, H_2D_list ,S_2D_list  = calc_step_constants(thermo_total_2D,keys,T,micromol,hr) 
    
    #k_f_HO, k_r_HO, G_HO_list, S_HO_list, H_HO_list = calc_step_constants(thermo_total_2D,keys,T,micromol,hr)
    
    
    '''
    ODE Solution Begins Here
    '''
    # =============================================================================
    # ODE params
    # =============================================================================
    
    #print('Number of reaction steps: ',len(TST_keys))
    #print('Number of gas species: ',len(gas_min_keys))
    #print('Number of surface species: ',len(local_min_keys))
    
    #The following parameters should be manually adjusted based on the ordering found within your
    #"Shomate.xlsx" spreadsheet.
    
    #The "shift" variable shifts the ordering for gas-phase flowrates.
    #Ensure that len(keys[-1]) is the number of stable species.-Not TSTs! (It's basically the number of species in the K column within "Shomate.xlsx")
    
    #Because the order for gas-phase begins at 0 (i.e "0g"), and the "F" vector contains both, ads and gas -phase species (in that specific order),
    #the "shift" variable is used to shift the ordering for F. For example, if one of the loops reaches "0g", its place in the "F" vector
    #is actually F[0+shift]. Or, once it reaches "1g", its place in the "F" vector is F[1+shift].
    shift = (len(keys[-1]))
    
    #Initial conditions, see "MKM_dir_list.xlsx" for numbering scheme.
    #For reference, F0[0:24] are ads species, F0[25:32] are gas species.
    F0 = np.zeros((len(local_min_keys))+len(gas_min_keys)+1)#Plus 1 because inert. Do not change this variable!
    F0[0] = 1/(1+(k_f_HO[0]/k_r_HO[0])) #Initial condition for O17 based on equilibrium constant.
    F0[1] = F0[0]*(k_f_HO[0]/k_r_HO[0]) #Initial condition for O16 see above.
    F0[-1] = (1/R/T)*v/micromol #Initial condition for inert efluent to keep reactor at 1 [atm]. Flowrate is [micromol/hr]
    if F0[-1] < 0:
        print("Warning, Inlet flowrate of inert is negative. Lower reactant concentration.")
        sys.exit()
    
    #Concentration of influent.
    #C_A_in = (1.0E-9)/micromol #[micromol/L]. Otherwise [mol/L]
    C_A_in = (0.004)/micromol #[micromol/L]. Otherwise [mol/L]
    C_B_trans_in = ((C_A_in/4)/(1+(k_f_HO[21]/k_r_HO[21]))) #[micromol/L]. Otherwise [mol/L]
    C_B_cis_in = C_B_trans_in*(k_f_HO[21]/k_r_HO[21]) #[micromol/L]. Otherwise [mol/L]
    
    #Flowrates into reactor.
    F_in = np.zeros((len(local_min_keys))+len(gas_min_keys)+1)#Plus 1 because inert. Do not change this variable!
    F_in[18] = v*C_A_in #gas-phase ethylene influent [micromol/hr]. Otherwise [mol/s]
    F_in[19] = v*C_B_trans_in #gas-phase trans-isoprene influent [micromol/hr]. Otherwise [mol/s]
    F_in[25] = v*C_B_cis_in #gas-phase cis-isoprene influent [micromol/hr]. Otherwise [mol/s]
    F_in[26] = (1-((C_A_in+C_B_trans_in+C_B_cis_in)*micromol*R*T))/R/T*v/micromol #Ensure inert flowrate keeps feed at 1 [atm] partial pressure.
    
    # ODE solver parameters
    abserr = 1.0e-30
    relerr = 1.0e-8
    
    stoptime = 500000#[hr]
    if micromol == 1:
        stoptime = 2000000000 #[hr], otherwise [s]
    
    
    #t1 = time.time()
    sol = solve_ivp(vectorfield, [0,stoptime], F0, method='BDF', args=([k_f_HO,k_r_HO],keys),rtol=relerr,atol=abserr,t_eval=np.linspace(0,stoptime,stoptime,endpoint=True))
    #t2 = time.time()
    #print('\nODE solver duration: %s' % (t2-t1))
    
    #Collect time
    t = sol.t
    
    #Collect solution & transpose.
    wsol = np.transpose(sol.y.copy())
    
    #Sanity check. Make sure site balance is preserved, no negative concetrations and pressure ~ 1[atm].
    for i,F_sol in enumerate(wsol[:,0]):
        sanity_check(wsol[i,:], shift, T, v, micromol,abserr,relerr,t[i])
    
    #Check if steady-state was approximately reached.
    SS_tol = 0.00001
    for i,species in enumerate(wsol[-1,:]):
        if np.absolute(species-wsol[-10,i]) > SS_tol:
            print("Warning: SS may not have been reached, try larger stoptime & plot transient flowrates.")
            sys.exit()
            
    print('Perturbed Vector')
    print(h_)#Print perturbation.
    print(h_intermediates)
    
    '''
    Analysis of Results
    '''
    
    #Conversion ethylene and total isoprene respectively:
    x_A = (F_in[18]-wsol[-1,18])/F_in[18]
    x_B = ((F_in[19]+F_in[25])-(wsol[-1,19]+wsol[-1,25]))/(F_in[19]+F_in[25])
    
    #Ensure conversion is ~1%
    if x_B > 0.01 or x_A > 0.01:
        print("Warning: Conversion is larger than ~1%: ",x_B,x_A)
    
    
    '''
    Print Results
    '''
    

    print('Reaction Rate [micro-mol/hr]')
    print('#','step','rate')
    for i,r in enumerate(rates[-1]):
        reactant_indices, product_indices, gas_phase_reaction = parse_key(keys[1][i])
        if gas_phase_reaction:
            print(i,keys[1][i],rates[-1][i]*V_g)
        else:
            print(i,keys[1][i],rates[-1][i]*W*Sigma)#*micromol/hr*R*T)

    #return rates[-1][11]#Change this to change "overall rxn rate".

if __name__ == "__main__":
    #Perform central finite difference to calculate Degree of rate control.
    #Let "h" be the perturbation.- Must be list of size related to steps.
    h_ = np.zeros((32))#There are 32 elementary steps in our reaction network.
    h_intermediates = np.zeros((18))#There are 18 intermediates.
    
    
    perturb = np.random.uniform(-0.05,0.05,5)
    #perturb = 0.001#Degree of perturbation.
    shift = 11 #DA product formation steps.
    #shift = 27 #DA product ads/des steps.
    #Iterate over ads elementary steps only!.
    for i in range(len(h_[shift:16])): #shift:16 for DA formation steps & shift:32 for product ads/des steps.
        h_[i+shift] = perturb[i]
        #r_foreward = CSTR_MKM(h_)#"overall rate" was chosen to be production of C7.
        #r_backward = CSTR_MKM(-h_)
        
        #X_RC = -R_J/1000*368.15*(np.log(r_foreward)-np.log(r_backward))/(perturb)
        
        #print(i+shift,X_RC)
        
       
    perturb = np.random.uniform(-0.05,0.05,5)
    shift = 27 #DA product formation steps.
    #shift = 27 #DA product ads/des steps.
    #Iterate over ads elementary steps only!.
    for i in range(len(h_[shift:32])): #shift:16 for DA formation steps & shift:32 for product ads/des steps.
        h_[i+shift] = perturb[i]    
    

    perturb = np.random.uniform(-0.05,0.05,5)
    shift = 13 #DA product intermediates.
    #Iterate over cycloadduct intermediates only!.
    for i in range(len(h_intermediates[shift:])): #shift: for DA cycloadducts as adsorbed intermedietes.
        h_intermediates[i+shift] = perturb[i]    
        
    h_intermediates[8] = np.random.uniform(-0.05,0.05,1)#For perturbing adsorbed trans-isoprene @ O17.
    h_intermediates[2] = np.random.uniform(-0.05,0.05,1)#For perturbing adsorbed ethene @ O17.

    CSTR_MKM(h_,h_intermediates)