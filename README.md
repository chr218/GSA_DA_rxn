# GSA_DA_rxn
This repo illustrates how to perform an ANOVA-based Global Sensitivity Analysis on the [4+2] Diels-Alder reactions between ethene and isoprene within H-ZSM5.

The two folders contain: 
1) "GSA_DA_rxn/MKM/" The microkinetic model of the [4+2] Diels-Alder reactions between ethene and isoprene within H-ZSM5, modified to generate predicted reaction rates under a 
perturbation of the TST barriers and intermediate free energies (and thereby the respective rate constants) from their nominal condition. 

2) "GSA_DA_rxn/Global_Sensitivity_Analysis/" The ANOVA-based Global Sensitivty Analysis (GSA), which fits the perturbed rate data from the MKM against the perturbation vector itself by applying 
polynomial chaos expansion (PCE). The polynomials from PCE are indepdent, which allows one to calculate the Sobol indices[1] (sensitivities) of each reaction
step & intermediate.

This approach is based on the work outlined by Tian and Rangarajan [2].

1. Sobol, I. M. Global sensitivity indices for nonlinear mathematical models and their Monte Carlo estimates. Math. Comput. Simulat. 2001, 55, 271– 280
2. Tian, H.; Rangarajan, S. Computing a Global Degree of Rate Control for Catalytic Systems. ACS Catalysis 2020, 10, 13535–13542. 


