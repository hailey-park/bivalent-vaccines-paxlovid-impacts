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
simulated_vax_params <- readRDS('data/simulated-ve-parameters.RDS')
bivalent_booster_data <- read.csv("outcome data clean/bivalent_booster_doses_clean_new_v2.csv")[,-1]
simulated_vax_waning_params <- readRDS('data/simulated-ve-waning-parameters.RDS')
bivalent_coverage <- read.csv("outcome data clean/bivalent_coverage_age_vax.csv")[,-1]


#Number of simulations
n <- 1000


####################################################
#Paxlovid/Bivalent Simulations

#Create df to store results
pax_hosp_results <- data.frame(matrix(nrow = 1000, ncol = 2*4))
colnames(pax_hosp_results) <- c(sprintf("strat%d_perc", 3:6),
                                 sprintf("strat%d_total", 3:6))
pax_death_results <- data.frame(matrix(nrow = 1000, ncol = 2*4))
colnames(pax_death_results) <- c(sprintf("strat%d_perc", 3:6),
                                sprintf("strat%d_total", 3:6))

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
  
  #create ve parameter df
  ve <- simulated_vax_params[i,]
  ve_estimates <- data.frame(baseline_vacc = c("Unvaccinated", "Primary Series", "Boosted (1 dose)","Boosted (2 doses)"),
                             ve_hosps = c(ve$unvax_hosp, ve$prim_hosp, ve$boost_hosp, ve$boost2_hosp),
                             ve_deaths = c(ve$unvax_death, ve$prim_death, ve$boost_death, ve$boost2_death))
  
  #merge together
  predictions_pax_df <- merge(merge(merge(covid_hosp 
                                          %>% group_by(age_group, boost_2_vax_status_match, vacc_status) %>% 
                                            summarise(predicted_cases = sum(predicted_cases),
                                                      predicted_hosp = sum(predicted_hosp)),
                                          pe_estimates, by.x = "vacc_status", by.y = "baseline_vacc", all.x = TRUE),
                                    ve_estimates, by.x = "boost_2_vax_status_match", by.y = "baseline_vacc", all.x = TRUE),
                              bivalent_coverage, by = c("boost_2_vax_status_match", "age_group"), all.x = TRUE) %>%
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
  
  #Strategy 3: Giving Paxlovid to those 18+ years (unvaccinated and comorbidities)
  strat3 <- predictions_pax_df %>% mutate(averted_hosp_diff_vax = case_when((age_group %in% c("18-49 years","50-64 years", "65-74 years", "75-84 years", "85+ years")) ~ predicted_hosp * ve_hosps * (0.7 - bivalent_coverage),
                                                                             TRUE ~ 0), 
                                          averted_hosp_diff_pax = case_when((age_group %in% c("18-49 years","50-64 years", "65-74 years", "75-84 years", "85+ years")) ~ (predicted_hosp - averted_hosp_diff_vax) * pe_hosps * (max_pax_uptake - baseline_pax_uptake),
                                                                             TRUE ~ 0),
                                          total_predicted = predicted_hosp * (1 - (baseline_pax_uptake * pe_hosps) - (bivalent_coverage * ve_hosps)))
  
  
  strat3$strategy <- rep("strat3", nrow(strat3))
  
  #Strategy 4: Giving Paxlovid to those 50+ years (unvaccinated and comorbidities)
  strat4 <- predictions_pax_df %>% mutate(averted_hosp_diff_vax = case_when((age_group %in% c("50-64 years", "65-74 years", "75-84 years", "85+ years")) ~ predicted_hosp * ve_hosps * (0.7 - bivalent_coverage),
                                                                             TRUE ~ 0),
                                          averted_hosp_diff_pax = case_when((age_group %in% c("50-64 years", "65-74 years", "75-84 years", "85+ years")) ~ (predicted_hosp - averted_hosp_diff_vax) * pe_hosps * (max_pax_uptake - baseline_pax_uptake),
                                                                             TRUE ~ 0),
                                          total_predicted = predicted_hosp * (1 - (baseline_pax_uptake * pe_hosps) - (bivalent_coverage * ve_hosps)))  
  
  strat4$strategy <- rep("strat4", nrow(strat4))
  
  
  #Strategy 5: Giving Paxlovid to the those 65+ years
  strat5 <- predictions_pax_df %>% mutate(averted_hosp_diff_vax = case_when((age_group %in% c("65-74 years", "75-84 years", "85+ years")) ~ predicted_hosp * ve_hosps * (0.7 - bivalent_coverage),
                                                                             TRUE ~ 0),
                                          averted_hosp_diff_pax = ifelse((age_group %in% c("65-74 years", "75-84 years", "85+ years")),
                                                                          (predicted_hosp - averted_hosp_diff_vax) * pe_hosps * (max_pax_uptake - baseline_pax_uptake), 0),
                                          total_predicted = predicted_hosp * (1 - (baseline_pax_uptake * pe_hosps) - (bivalent_coverage * ve_hosps)))
  
  strat5$strategy <- rep("strat5", nrow(strat5))
  
  #Strategy 6: Giving Paxlovid to the those 75+ years
  strat6 <- predictions_pax_df %>% mutate(averted_hosp_diff_vax = case_when((age_group %in% c("75-84 years", "85+ years")) ~ predicted_hosp * ve_hosps * (0.7 - bivalent_coverage),
                                                                             TRUE ~ 0),
                                          averted_hosp_diff_pax = ifelse((age_group %in% c("75-84 years", "85+ years")),
                                                                          (predicted_hosp - averted_hosp_diff_vax) * pe_hosps * (max_pax_uptake - baseline_pax_uptake), 0),
                                          total_predicted = predicted_hosp * (1 - (baseline_pax_uptake * pe_hosps) - (bivalent_coverage * ve_hosps)))
  
  strat6$strategy <- rep("strat6", nrow(strat6))  
  
  combined_strat <- bind_rows(strat3, strat4, strat5, strat6)
  
  
  #######################################
  combined_strat_outcomes <- combined_strat %>% group_by(strategy) %>% summarise(averted_outcomes_diff_vax = ceiling(sum(averted_hosp_diff_vax)),
                                                                                 averted_outcomes_diff_pax = ceiling(sum(averted_hosp_diff_pax)),
                                                                                 total_hosp = ceiling(sum(total_predicted))) %>%
    mutate(averted_outcomes_cascade_pax = ceiling(averted_outcomes_diff_pax * 0.8 * 0.895),
           total_averted_outcomes = averted_outcomes_cascade_pax + averted_outcomes_diff_vax,
           perc_averted = round(total_averted_outcomes/ total_hosp * 100, 1),
           strategy = factor(strategy, levels = c("strat3", "strat4", "strat5", "strat6")))  %>%
    arrange(strategy)
  
  pax_hosp_results[i,] <- c(combined_strat_outcomes$perc_averted, combined_strat_outcomes$total_averted_outcomes)
}

saveRDS(pax_hosp_results, "data/results-combined-hosp.RDS")

#################################################################
#Evaluate the 95% UIs for the results
results_combo_death <- readRDS("data/results-combined-death.RDS")
results_combo_hosp <- readRDS("data/results-combined-hosp.RDS")




inspection <- results_combo_hosp %>% apply(2, quantile, probs=c(0.025, 0.975))





