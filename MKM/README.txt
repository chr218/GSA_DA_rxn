Author: Christopher Rzepa

These files contain the MKM for the [4+2] Diels-Alder reaction
between ethene and isoprene within H-ZSM5.

File "CSTR_MKM_Lite_ordered.py" is the microkinetic model. The
thermodynamic parameters (Gibbs free energies of TST and intermediates)
are calculated within the "data_config.py" file. The Shomate parameters
for these species are held within "Shomate_ordered.xlsx".

To perform an ANOVA-baed Global Sensitivity Analysis, one should
run this MKM ~5,000 times. An example of the output is provided in
.csv format for perturbations of 0.5 [kJ/mol] as
"output_0.5.csv.zip" etc.
