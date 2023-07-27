###################################################################################################
#Title: Bivalent Booster Regression Analysis (with 95% UI)
#Author: Hailey Park
#Date: January 10, 2023
###################################################################################################

rm(list=ls())

setwd("/mnt/projects/covid_partners/ucsf_lo/Bivalent Booster Paxlovid Impact Analysis/final-code")


#Loading in libraries
library(readr)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(tidyr)
library(tibble)
library(reshape2)
library(lubridate)
library(scales)
library(data.table)
library(ciTools)

#Read in dataset
base_df <- readRDS("data/base-dataframe.RDS") 
simulated_cases <- readRDS('data/simulated-cases.RDS')
simulated_hosp <- readRDS('data/simulated-hosp.RDS')
simulated_death <- readRDS('data/simulated-death.RDS')
simulated_pax_params <- readRDS('data/simulated-pax-parameters.RDS')

#Number of simulations
n <- 1000

####################################################
#Paxlovid Simulations

#Create df to store results
pax_hosp_results <- data.frame(matrix(nrow = 1000, ncol = 3*2))
colnames(pax_hosp_results) <- c(sprintf("strat%d_nnt", 1:2),
                                 sprintf("strat%d_perc", 1:2),
                                 sprintf("strat%d_total", 1:2))
pax_death_results <- data.frame(matrix(nrow = 1000, ncol = 3*2))
colnames(pax_death_results) <- c(sprintf("strat%d_nnt", 1:2),
                                sprintf("strat%d_perc", 1:2),
                                sprintf("strat%d_total", 1:2))

#Run Monte Carlo simulation (N = 1000)
for (i in c(1:n)) {
  
  print(i)

  #create df with predictions
  covid_hosp <- base_df %>%  mutate(vacc_status = ifelse(boost_2_vax_status_match %in% c("Primary Series", "Boosted (1 dose)", "Boosted (2 doses)"),
                                                                                                      "Vaccinated",
                                                                                                      "Unvaccinated"))
  covid_hosp$predicted_hosp <- simulated_hosp[i,]
  covid_hosp$predicted_cases <- simulated_cases[i,]
  
  #create pe parameter df
  pe <- simulated_pax_params[i,]
  pe_estimates <- data.frame(baseline_vacc = c("Unvaccinated", "Vaccinated"),
                             pe_hosps = c(pe$unvax_hosp, pe$vax_hosp),
                             pe_deaths = c(pe$unvax_death, pe$vax_death))
  
  #merge together
  predictions_pax_df <- merge(covid_hosp 
                              %>% group_by(age_group, boost_2_vax_status_match, vacc_status) %>% 
                                summarise(predicted_cases = sum(predicted_cases),
                                          predicted_hosp = sum(predicted_hosp)),
                              pe_estimates, by.x = "vacc_status", by.y = "baseline_vacc", all.x = TRUE) %>%
    mutate(baseline_pax_uptake = case_when(age_group == "50-64 years" ~ 0.237,
                                           age_group == "65-74 years" ~ 0.341,
                                           age_group == "75-84 years" ~ 0.332,
                                           age_group == "85+ years" ~ 0.322,
                                           TRUE ~ 0),
           max_pax_uptake = case_when(age_group == "18-49 years" ~ 0.02,
                                      age_group == "50-64 years" ~ 0.33,
                                      age_group %in% c("65-74 years", "75-84 years", "85+ years") ~ 0.7,
                                      TRUE ~ 0))
  
  #############
  #simulate strategies
  #Strategy 1: Giving Paxlovid to those 18+ years (everyone)
    #Assumption: setting maximum paxlovid coverage to be 70% for everyone
  strat1 <- predictions_pax_df %>% mutate(averted_hosp_diff = case_when((age_group %in% c("18-49 years", "50-64 years", "65-74 years", "75-84 years", "85+ years"))  ~ predicted_hosp * pe_hosps * (0.7 - baseline_pax_uptake),
                                                                         TRUE ~ 0),
                                          target_group_cases_diff = case_when((age_group %in% c("18-49 years", "50-64 years", "65-74 years", "75-84 years", "85+ years")) ~ predicted_cases * (0.7 - baseline_pax_uptake),
                                                                              TRUE ~ 0),
                                          total_predicted = predicted_hosp * (1 - (baseline_pax_uptake * pe_hosps)))
  strat1$strategy <- rep("strat1", nrow(strat1))
  
  #Strategy 2: Giving Paxlovid to those 50+ years (everyone)
    #Assumption: setting maximum paxlovid coverage to be 70% for everyone
  strat2 <- predictions_pax_df %>% mutate(averted_hosp_diff = case_when((age_group %in% c("50-64 years", "65-74 years", "75-84 years", "85+ years"))  ~ predicted_hosp * pe_hosps * (0.7 - baseline_pax_uptake),
                                                                         TRUE ~ 0),
                                          target_group_cases_diff = case_when((age_group %in% c("50-64 years", "65-74 years", "75-84 years", "85+ years")) ~ predicted_cases * (0.7 - baseline_pax_uptake),
                                                                              TRUE ~ 0),
                                          total_predicted = predicted_hosp * (1 - (baseline_pax_uptake * pe_hosps)))
  strat2$strategy <- rep("strat2", nrow(strat2))
  
  combined_strat <- bind_rows(strat1, strat2)
  
  
  #######################################
  combined_strat_outcomes <- combined_strat %>% group_by(strategy) %>% summarise(averted_outcomes_diff = ceiling(sum(averted_hosp_diff)),
                                                                                 total_hosp = ceiling(sum(total_predicted)),
                                                                                 target_group_cases_diff = ceiling(sum(target_group_cases_diff))) %>%
    mutate(averted_outcomes_cascade = ceiling(averted_outcomes_diff * 0.8 * 0.895),
           treated_cascade = ceiling(target_group_cases_diff * 0.8 * 0.895),
           NNT = ceiling(treated_cascade/averted_outcomes_cascade),
           perc_averted = round(averted_outcomes_cascade/ total_hosp * 100, 1),
           strategy = factor(strategy, levels = c("strat1", "strat2"))) %>%
    arrange(strategy)
  
  pax_hosp_results[i,] <- c(combined_strat_outcomes$NNT, combined_strat_outcomes$perc_averted, combined_strat_outcomes$averted_outcomes_cascade)
}

saveRDS(pax_hosp_results, "data/results-sensitivity-expanded eligibility-pax-hosp.RDS")

#################################################################
#Evaluate the 95% UIs for the results
results_pax_death <- readRDS("data/results-sensitivity-expanded eligibility-pax-death.RDS")
results_pax_hosp <- readRDS("data/results-sensitivity-expanded eligibility-pax-hosp.RDS")



inspection <- results_pax_hosp %>% apply(2, quantile, probs=c(0.025, 0.975))





