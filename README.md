# Predicting the Public Health Impact of Bivalent Vaccines and Nirmatrelvir-Ritonavir Against COVID-19
This repository contains analytic code for modeling the impact of increased uptake of bivalent COVID-19 vaccines and nirmatrelvir-ritonavir (Paxlovid) during acute illness in different risk groups.

Data used for this analysis is not publicly available. Data requests may be made to the California Department of Public Health. 


## Structure
* `1-data cleaning`: contains all code needed for initial cleaning of data
* `2-analysis`:
  * `1-waning model`: contains code for constructing the waning protection curves by risk group, and waning curves used for sensitivity analyses
  * `2-model calibration`: contains code for constructing hypothetical cohorts of 1 million individuals for each risk group and calibrating age-specific estimates of severe COVID-19 risk 
  * `3-simulation main`: contains code for running different frequencies of bivalent booster vaccination (no booster, one-time booster, annual booster, semiannual booster) over a two year time horizon
  * `4-simulation five years`: contains code for running different frequencies of bivalent booster vaccination over a five year time horizon
  *  `5-simulation delayed vaccine admin`: contains code for running different frequencies of bivalent booster vaccination with delayed vaccine rollout (6-months)
  *  `6-simulation lower ve after 1st dose`: contains code for running different frequencies of bivalent booster vaccination with assumption of lower vaccine effectiveness after the 1st bivalent booster dose
  *  `7-model validation`: contains code for running a model validation over a 3-month period 
* `figures tables`: contains figure and table outputs (jpg and csv)
