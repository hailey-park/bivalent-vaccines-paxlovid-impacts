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
bivalent_booster_data <- read.csv("outcome data clean/bivalent_booster_doses_clean_sensitivity_reverse_uptake.csv")[,-1]
simulated_vax_waning_params <- readRDS('data/simulated-ve-waning-parameters.RDS')
bivalent_coverage <- read.csv("outcome data clean/bivalent_coverage_age_vax.csv")[,-1]


#Number of simulations
n <- 1000

####################################################
#Bivalent Vaccine Simulations

#Create df to store results
vax_hosp_results <- data.frame(matrix(nrow = 1000, ncol = 3*7))
colnames(vax_hosp_results) <- c(sprintf("strat%d_nnt", 1:7),
                                 sprintf("strat%d_perc", 1:7),
                                 sprintf("strat%d_total", 1:7))
vax_death_results <- data.frame(matrix(nrow = 1000, ncol = 3*7))
colnames(vax_death_results) <- c(sprintf("strat%d_nnt", 1:7),
                                sprintf("strat%d_perc", 1:7),
                                sprintf("strat%d_total", 1:7))

#Run Monte Carlo simulation (N = 1000)
for (i in c(1:n)) {
  print(i)
  #create df with predictions
  covid_df <- base_df 
  covid_df$predicted_outcomes <- simulated_hosp[i,]
  
  #create ve parameter df
  ve <- simulated_vax_params[i,]
  ve_estimates <- data.frame(baseline_vacc = c("Unvaccinated", "Primary Series", "Boosted (1 dose)","Boosted (2 doses)"),
                             ve_hosps = c(ve$unvax_hosp, ve$prim_hosp, ve$boost_hosp, ve$boost2_hosp),
                             ve_deaths = c(ve$unvax_death, ve$prim_death, ve$boost_death, ve$boost2_death))
  
  #merge together
  predictions_df <- merge(merge(covid_df %>% group_by(age_group, boost_2_vax_status_match) %>% 
                            summarise(predicted_outcomes = sum(predicted_outcomes)),
                          ve_estimates, by.x = "boost_2_vax_status_match", by.y = "baseline_vacc", all.x = TRUE),
                          bivalent_coverage, by = c("boost_2_vax_status_match", "age_group"), all.x = TRUE) %>%
    mutate(predicted_outcomes = round(predicted_outcomes),
           max_uptake = case_when(age_group %in% c("0-17 years", "18-49 years", "50-64 years") ~ 0.7,
                                  age_group %in% c("65-74 years", "75-84 years", "85+ years") ~ 0.5,
                                  TRUE ~ 0),
           max_uptake = ifelse(bivalent_coverage > max_uptake, bivalent_coverage, max_uptake))
  
  #############
  #simulate strategies
  strat1 <- predictions_df %>% mutate(averted_diff = ifelse((boost_2_vax_status_match %in% c("Unvaccinated", "Primary Series", "Boosted (1 dose)", "Boosted (2 doses)")),
                                                            predicted_outcomes * ve_hosps * (max_uptake - bivalent_coverage), 0),
                                      total_predicted = predicted_outcomes * (1 - (bivalent_coverage * ve_hosps)))
  
  strat1$strategy <- rep("strat1", nrow(strat1))
  
  strat2 <- predictions_df %>% mutate(averted_diff = ifelse((boost_2_vax_status_match %in% c("Primary Series", "Boosted (1 dose)", "Boosted (2 doses)")),
                                                            predicted_outcomes * ve_hosps * (max_uptake - bivalent_coverage), 0),
                                      total_predicted = predicted_outcomes * (1 - (bivalent_coverage * ve_hosps)))
  strat2$strategy <- rep("strat2", nrow(strat2))
  
  strat3 <- predictions_df %>% mutate(averted_diff = ifelse((boost_2_vax_status_match %in% c("Unvaccinated")),
                                                            predicted_outcomes * ve_hosps * (max_uptake - bivalent_coverage), 0),
                                      total_predicted = predicted_outcomes * (1 - (bivalent_coverage * ve_hosps)))
  strat3$strategy <- rep("strat3", nrow(strat3))
  
  strat4 <- predictions_df %>% mutate(averted_diff = ifelse((boost_2_vax_status_match %in% c("Primary Series")),
                                                            predicted_outcomes * ve_hosps * (max_uptake - bivalent_coverage), 0),
                                      total_predicted = predicted_outcomes * (1 - (bivalent_coverage * ve_hosps)))
  strat4$strategy <- rep("strat4", nrow(strat3))
  
  strat5 <- predictions_df %>% mutate(averted_diff = ifelse(((age_group %in% c( "75-84 years", "85+ years")) & (boost_2_vax_status_match != "Unvaccinated")),
                                                            predicted_outcomes * ve_hosps * (max_uptake - bivalent_coverage), 0),
                                      total_predicted = predicted_outcomes * (1 - (bivalent_coverage * ve_hosps))) 
  strat5$strategy <- rep("strat5", nrow(strat3))
  
  strat6 <- predictions_df %>% mutate(averted_diff = ifelse(((age_group %in% c( "65-74 years", "75-84 years", "85+ years")) & (boost_2_vax_status_match != "Unvaccinated")),
                                                            predicted_outcomes * ve_hosps * (max_uptake - bivalent_coverage), 0),
                                      total_predicted = predicted_outcomes * (1 - (bivalent_coverage * ve_hosps)))
  strat6$strategy <- rep("strat6", nrow(strat1))
  
  strat7 <- predictions_df %>% mutate(averted_diff = ifelse(((age_group %in% c("50-64 years", "65-74 years", "75-84 years", "85+ years")) & (boost_2_vax_status_match != "Unvaccinated")),
                                                            predicted_outcomes * ve_hosps * (max_uptake - bivalent_coverage), 0),
                                      total_predicted = predicted_outcomes * (1 - (bivalent_coverage * ve_hosps)))
  strat7$strategy <- rep("strat7", nrow(strat1))
  
  
  
  combined_strat <- bind_rows(strat1, strat2, strat3, strat4, strat5, strat6, strat7)
  
  
  combined_strat_outcomes <- combined_strat %>% group_by(strategy) %>% summarise(averted_outcomes_diff = round(sum(averted_diff)),
                                                                                 total_predicted_outcomes = round(sum(total_predicted))) 
  
  
  merged_data <- merge(combined_strat_outcomes, bivalent_booster_data, by = c("strategy"), all.X = TRUE) %>%
    mutate(ARR = averted_outcomes_diff/(total_doses),
           NNT = ceiling(1/ARR),
           perc_averted = round(averted_outcomes_diff/total_predicted_outcomes * 100, 1),
           strategy = factor(strategy, levels = c("strat1", "strat2", "strat3", "strat4", "strat5", "strat6", "strat7")))%>%
    select(-c(variable)) %>% arrange(strategy)
  
  vax_hosp_results[i,] <- c(merged_data$NNT, merged_data$perc_averted, merged_data$averted_outcomes_diff)
}

saveRDS(vax_hosp_results, "data/results-sensitivity-reverseUptake-vax-hosp.RDS")
####################################################
#Paxlovid Simulations

#Create df to store results
pax_hosp_results <- data.frame(matrix(nrow = 1000, ncol = 3*4))
colnames(pax_hosp_results) <- c(sprintf("strat%d_nnt", 3:6),
                                 sprintf("strat%d_perc", 3:6),
                                 sprintf("strat%d_total", 3:6))
pax_death_results <- data.frame(matrix(nrow = 1000, ncol = 3*4))
colnames(pax_death_results) <- c(sprintf("strat%d_nnt", 3:6),
                                sprintf("strat%d_perc", 3:6),
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
                                      age_group %in% c("65-74 years", "75-84 years", "85+ years") ~ 0.5,
                                      TRUE ~ 0))
  
  #############
  #simulate strategies
  
  #Strategy 3: Giving Paxlovid to those 18+ years (unvaccinated and comorbidities)
  strat3 <- predictions_pax_df %>% mutate(averted_hosp_diff = case_when((age_group %in% c("18-49 years","50-64 years", "65-74 years", "75-84 years", "85+ years")) ~ predicted_hosp * pe_hosps * (max_pax_uptake - baseline_pax_uptake),
                                                                         TRUE ~ 0),
                                          target_group_cases_diff = case_when((age_group %in% c("18-49 years","50-64 years", "65-74 years", "75-84 years", "85+ years")) ~ predicted_cases * (max_pax_uptake - baseline_pax_uptake),
                                                                              TRUE ~ 0),
                                          total_predicted = predicted_hosp * (1 - (baseline_pax_uptake * pe_hosps)))
  
  
  strat3$strategy <- rep("strat3", nrow(strat3))
  
  #Strategy 4: Giving Paxlovid to those 50+ years (unvaccinated and comorbidities)
  strat4 <- predictions_pax_df %>% mutate(averted_hosp_diff = case_when((age_group %in% c("50-64 years", "65-74 years", "75-84 years", "85+ years")) ~ predicted_hosp * pe_hosps * (max_pax_uptake - baseline_pax_uptake),
                                                                         TRUE ~ 0),
                                          target_group_cases_diff = case_when((age_group %in% c("50-64 years", "65-74 years", "75-84 years", "85+ years")) ~ predicted_cases * (max_pax_uptake - baseline_pax_uptake),
                                                                              TRUE ~ 0),
                                          total_predicted = predicted_hosp * (1 - (baseline_pax_uptake * pe_hosps)))  
  
  strat4$strategy <- rep("strat4", nrow(strat4))
  
  
  #Strategy 5: Giving Paxlovid to the those 65+ years
  strat5 <- predictions_pax_df %>% mutate(averted_hosp_diff = ifelse((age_group %in% c("65-74 years", "75-84 years", "85+ years")),
                                                                      predicted_hosp * pe_hosps * (max_pax_uptake - baseline_pax_uptake), 0),
                                          target_group_cases_diff = ifelse((age_group %in% c("65-74 years", "75-84 years", "85+ years")),
                                                                           predicted_cases * (max_pax_uptake - baseline_pax_uptake), 0),
                                          total_predicted = predicted_hosp * (1 - (baseline_pax_uptake * pe_hosps)))
  
  strat5$strategy <- rep("strat5", nrow(strat5))
  
  #Strategy 6: Giving Paxlovid to the those 75+ years
  strat6 <- predictions_pax_df %>% mutate(averted_hosp_diff = ifelse((age_group %in% c("75-84 years", "85+ years")),
                                                                      predicted_hosp * pe_hosps * (max_pax_uptake - baseline_pax_uptake), 0),
                                          target_group_cases_diff = ifelse((age_group %in% c("75-84 years", "85+ years")),
                                                                           predicted_cases * (max_pax_uptake - baseline_pax_uptake), 0),
                                          total_predicted = predicted_hosp * (1 - (baseline_pax_uptake * pe_hosps)))
  
  strat6$strategy <- rep("strat6", nrow(strat6))  
  

  
  
  combined_strat <- bind_rows(strat3, strat4, strat5, strat6)
  
  
  #######################################
  combined_strat_outcomes <- combined_strat %>% group_by(strategy) %>% summarise(averted_outcomes_diff = ceiling(sum(averted_hosp_diff)),
                                                                                 total_hosp = ceiling(sum(total_predicted)),
                                                                                 target_group_cases_diff = ceiling(sum(target_group_cases_diff))) %>%
    mutate(averted_outcomes_cascade = ceiling(averted_outcomes_diff * 0.8 * 0.895),
           treated_cascade = ceiling(target_group_cases_diff * 0.8 * 0.895),
           NNT = ceiling(treated_cascade/averted_outcomes_cascade),
           perc_averted = round(averted_outcomes_cascade/ total_hosp * 100, 1),
           strategy = factor(strategy, levels = c("strat3", "strat4", "strat5", "strat6"))) %>%
    arrange(strategy)
  
  pax_hosp_results[i,] <- c(combined_strat_outcomes$NNT, combined_strat_outcomes$perc_averted, combined_strat_outcomes$averted_outcomes_cascade)
}

saveRDS(pax_hosp_results, "data/results-sensitivity-reverseUptake-pax-hosp.RDS")


#################################################################
#Waning Vaccine Simulations - Main Analysis for Cases

#Create df to store results
vax_cases_results <- data.frame(matrix(nrow = 1000, ncol = 3*7))
colnames(vax_cases_results) <- c(sprintf("strat%d_nnt", 1:7),
                                 sprintf("strat%d_perc", 1:7),
                                 sprintf("strat%d_total", 1:7))

#Run Monte Carlo simulation (N = 1000)
for (i in c(1:n)) {
  print(i)
  
  #create df with predictions
  covid_cases <- base_df 
  covid_cases$predicted_cases <- simulated_cases[i,]

  #create ve parameter df
  ve <- simulated_vax_waning_params[i,]

  #Add waning ve to df
  predictions_df <- merge(covid_cases %>% mutate(ve_cases = case_when((weeks_since_july2022 %in% c(0:4) & boost_2_vax_status_match == 'Unvaccinated') ~ ve$unvax_0_to_4,
                                                                (weeks_since_july2022 %in% c(5:26) & boost_2_vax_status_match == 'Unvaccinated') ~ ve$unvax_5_plus,
                                                                (weeks_since_july2022 %in% c(0:4) & boost_2_vax_status_match == 'Primary Series') ~ ve$prim_0_to_4,
                                                                (weeks_since_july2022 %in% c(5:9) & boost_2_vax_status_match == 'Primary Series') ~ ve$prim_5_to_9,
                                                                (weeks_since_july2022 %in% c(10:13) & boost_2_vax_status_match == 'Primary Series') ~ ve$prim_10_to_13,
                                                                (weeks_since_july2022 %in% c(14:17) & boost_2_vax_status_match == 'Primary Series') ~ ve$prim_14_to_17,
                                                                (weeks_since_july2022 %in% c(18:26) & boost_2_vax_status_match == 'Primary Series') ~ ve$prim_18_plus,
                                                                (weeks_since_july2022 %in% c(0:5) & boost_2_vax_status_match == 'Boosted (1 dose)') ~ ve$boost_0_to_5,
                                                                (weeks_since_july2022 %in% c(6:15) & boost_2_vax_status_match == 'Boosted (1 dose)') ~ ve$boost_6_to_15,
                                                                (weeks_since_july2022 %in% c(16:26) & boost_2_vax_status_match == 'Boosted (1 dose)') ~ ve$boost_16_plus,
                                                                (weeks_since_july2022 %in% c(0:5) & boost_2_vax_status_match == 'Boosted (2 doses)') ~ ve$boost_0_to_5,
                                                                (weeks_since_july2022 %in% c(6:15) & boost_2_vax_status_match == 'Boosted (2 doses)') ~ ve$boost_6_to_15,
                                                                (weeks_since_july2022 %in% c(16:26) & boost_2_vax_status_match == 'Boosted (2 doses)') ~ ve$boost_16_plus,
                                                                TRUE ~ 0)),
                          bivalent_coverage, by = c("boost_2_vax_status_match", "age_group"), all.x = TRUE) %>%
    mutate(predicted_cases = round(predicted_cases),
           max_uptake = case_when(age_group %in% c("0-17 years", "18-49 years", "50-64 years") ~ 0.7,
                                  age_group %in% c("65-74 years", "75-84 years", "85+ years") ~ 0.5,
                                  TRUE ~ 0),
           max_uptake = ifelse(bivalent_coverage > max_uptake, bivalent_coverage, max_uptake))
  
  #############
  #simulate strategies
  strat1 <- predictions_df %>% mutate(averted_diff = ifelse((boost_2_vax_status_match %in% c("Unvaccinated", "Primary Series", "Boosted (1 dose)", "Boosted (2 doses)")),
                                                            predicted_cases * ve_cases * (max_uptake - bivalent_coverage), 0),
                                      total_predicted = predicted_cases * (1 - (bivalent_coverage * ve_cases)))
  
  strat1$strategy <- rep("strat1", nrow(strat1))
  
  strat2 <- predictions_df %>% mutate(averted_diff = ifelse((boost_2_vax_status_match %in% c("Primary Series", "Boosted (1 dose)", "Boosted (2 doses)")),
                                                            predicted_cases * ve_cases * (max_uptake - bivalent_coverage), 0),
                                      total_predicted = predicted_cases * (1 - (bivalent_coverage * ve_cases)))
  strat2$strategy <- rep("strat2", nrow(strat2))
  
  strat3 <- predictions_df %>% mutate(averted_diff = ifelse((boost_2_vax_status_match %in% c("Unvaccinated")),
                                                            predicted_cases * ve_cases * (max_uptake - bivalent_coverage), 0),
                                      total_predicted = predicted_cases * (1 - (bivalent_coverage * ve_cases)))
  strat3$strategy <- rep("strat3", nrow(strat3))
  
  strat4 <- predictions_df %>% mutate(averted_diff = ifelse((boost_2_vax_status_match %in% c("Primary Series")),
                                                            predicted_cases * ve_cases * (max_uptake - bivalent_coverage), 0),
                                      total_predicted = predicted_cases * (1 - (bivalent_coverage * ve_cases)))
  strat4$strategy <- rep("strat4", nrow(strat3))
  
  strat5 <- predictions_df %>% mutate(averted_diff = ifelse(((age_group %in% c( "75-84 years", "85+ years")) & (boost_2_vax_status_match != "Unvaccinated")),
                                                            predicted_cases * ve_cases * (max_uptake - bivalent_coverage), 0),
                                      total_predicted = predicted_cases * (1 - (bivalent_coverage * ve_cases))) 
  strat5$strategy <- rep("strat5", nrow(strat3))
  
  strat6 <- predictions_df %>% mutate(averted_diff = ifelse(((age_group %in% c( "65-74 years", "75-84 years", "85+ years")) & (boost_2_vax_status_match != "Unvaccinated")),
                                                            predicted_cases * ve_cases * (max_uptake - bivalent_coverage), 0),
                                      total_predicted = predicted_cases * (1 - (bivalent_coverage * ve_cases)))
  strat6$strategy <- rep("strat6", nrow(strat1))
  
  strat7 <- predictions_df %>% mutate(averted_diff = ifelse(((age_group %in% c("50-64 years", "65-74 years", "75-84 years", "85+ years")) & (boost_2_vax_status_match != "Unvaccinated")),
                                                            predicted_cases * ve_cases * (max_uptake - bivalent_coverage), 0),
                                      total_predicted = predicted_cases * (1 - (bivalent_coverage * ve_cases)))
  strat7$strategy <- rep("strat7", nrow(strat1))
  
  
  
  combined_strat <- bind_rows(strat1, strat2, strat3, strat4, strat5, strat6, strat7)
  
  
  combined_strat_outcomes <- combined_strat %>% group_by(strategy) %>% summarise(averted_outcomes_diff = round(sum(averted_diff)),
                                                                                 total_predicted_outcomes = round(sum(total_predicted))) 
  
  
  merged_data <- merge(combined_strat_outcomes, bivalent_booster_data, by = c("strategy"), all.X = TRUE) %>%
    mutate(ARR = averted_outcomes_diff/(total_doses),
           NNT = ceiling(1/ARR),
           perc_averted = round(averted_outcomes_diff/total_predicted_outcomes * 100, 1),
           strategy = factor(strategy, levels = c("strat1", "strat2", "strat3", "strat4", "strat5", "strat6", "strat7")))%>%
    select(-c(variable)) %>% 
    arrange(strategy)
  
  
  vax_cases_results[i,] <- c(merged_data$NNT, merged_data$perc_averted, merged_data$averted_outcomes_diff)
}

saveRDS(vax_cases_results, "data/results-sensitivity-reverseUptake-vax-cases.RDS")

#################################################################
#Evaluate the 95% UIs for the results
results_pax_death <- readRDS("data/results-sensitivity-reverseUptake-pax-death.RDS")
results_pax_hosp <- readRDS("data/results-sensitivity-reverseUptake-pax-hosp.RDS")
results_vax_hosp <- readRDS("data/results-sensitivity-reverseUptake-vax-hosp.RDS")
results_vax_cases <- readRDS("data/results-sensitivity-reverseUptake-vax-cases.RDS")
results_vax_death <- readRDS("data/results-sensitivity-reverseUptake-vax-death.RDS")



inspection <- results_pax_hosp %>% apply(2, quantile, probs=c(0.025, 0.975))





