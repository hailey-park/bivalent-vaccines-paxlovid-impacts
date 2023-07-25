# Predicting the Public Health Impact of Bivalent Vaccines and Nirmatrelvir-Ritonavir Against COVID-19
This repository contains analytic code for modeling the impact of increased uptake of bivalent COVID-19 vaccines and nirmatrelvir-ritonavir (Paxlovid) during acute illness in different risk groups.

Data used for this analysis is not publicly available. Data requests may be made to the California Department of Public Health. 

This paper will be published in _Open Forum Infectious Diseases_.

## Structure
* `data cleaning`: contains all code needed for initial cleaning of data
* `analysis`:
  * `1-main analysis`: contains code for main regression analyses for bivalent vaccines and nirmatrelvir-ritonavir interventions
  * `2-cross validation`: contains code for running a cross-validation over a 3-month period 
  * `3-sensitivity-stratified risk`: contains code for regression analyses using risk groups stratified by age group and vaccination status
  * `4-sensitivity-expanded eligibility`: contains code for the nirmatrelvir-ritonavir regression analyses modeling expanded eligibility in the 18+ age group and 50+ age groups.
  *  `5-sensitivity-sustained background vaccination`: contains code for bivalent vaccine regression analyses modeling alternative status quo with increasing bivalent coverage over time.
  *  `6-sensitivity-2x case ascertainment`: contains code for running for bivalent vaccine regression analyses that accounts for low case ascertainment
  *  `7-sensitivity-perfect uptake`: contains code for running regression analyses assuming perfect (100%) uptake of interventions
  *  `8-sensitivity-uptake correlated with risk`: contains code for running regression analyses assuming uptake correlated with risk of COVID-19
  *  `9-sensitivity-uptake inversely correlated with risk`: contains code for running regression analyses assuming uptake inversely-correlated with risk of COVID-19
  *  `10-sensitivity-older age bivalent coverage`: contains code for running regression analyses modeling age-stratified baseline bivalent coverage in the oldest age groups.
  *  `11-sensitivity-combination paxlovid bivalent`: contains code for running regression analyses on the combined effect of both bivalent vaccine and nirmatrelvir-ritonavir interventions.
  *  `12-sensitivity-reduction in effectiveness`: contains code for comparison of NNT between high and low-risk groups on the reduction in effectiveness of bivalent vaccines and nirmatrelvir-ritonavir.
* `figures tables`: contains figure and table outputs (jpg and csv)
