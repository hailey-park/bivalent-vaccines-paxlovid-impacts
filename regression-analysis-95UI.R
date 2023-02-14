###################################################################################################
#Title: Bivalent Booster Regression Analysis (with 95% UI)
#Author: Hailey Park
#Date: January 10, 2023
###################################################################################################

rm(list=ls())

setwd("/mnt/projects/covid_partners/ucsf_lo/Bivalent Booster Paxlovid Impact Analysis")


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
simulated_pax_params <- readRDS('data/simulated-pax-parameters-real.RDS')
simulated_vax_params <- readRDS('data/simulated-ve-parameters.RDS')
bivalent_booster_data <- read.csv("outcome data clean/bivalent_booster_doses_clean_new.csv")[,-1]
simulated_vax_waning_params <- readRDS('data/simulated-ve-waning-parameters.RDS')
base_df <- readRDS("data/base-dataframe-crossval.RDS") 
crossval_simulated_cases <- readRDS('data/simulated-cases-crossval.RDS')
crossval_simulated_hosp <- readRDS('data/simulated-hosp-crossval.RDS')
crossval_simulated_death <- readRDS('data/simulated-death-crossval.RDS')

#Number of simulations
n <- 1000

####################################################
#Bivalent Vaccine Simulations

#Create df to store results
vax_cases_results <- data.frame(matrix(nrow = 1000, ncol = 3*16))
colnames(vax_cases_results) <- c(sprintf("strat%d_nnt", 1:16),
                                 sprintf("strat%d_perc", 1:16),
                                 sprintf("strat%d_total", 1:16))
vax_hosp_results <- data.frame(matrix(nrow = 1000, ncol = 3*16))
colnames(vax_hosp_results) <- c(sprintf("strat%d_nnt", 1:16),
                                sprintf("strat%d_perc", 1:16),
                                sprintf("strat%d_total", 1:16))
vax_death_results <- data.frame(matrix(nrow = 1000, ncol = 3*16))
colnames(vax_death_results) <- c(sprintf("strat%d_nnt", 1:16),
                                 sprintf("strat%d_perc", 1:16),
                                 sprintf("strat%d_total", 1:16))

#Run Monte Carlo simulation (N = 1000)
for (i in c(1:n)) {
  print(i)
  #create df with predictions
  covid_cases <- base_df 
  covid_cases$predicted_cases <- simulated_cases[i,]
  
  #create ve parameter df
  ve <- simulated_vax_params[i,]
  ve_estimates <- data.frame(baseline_vacc = c("Unvaccinated", "Primary Series", "Boosted (1 dose)","Boosted (2 doses)"),
                             ve_cases = c(ve$unvax_case, ve$prim_case, ve$boost_case, 0),
                             ve_hosps = c(ve$unvax_hosp, ve$prim_hosp, ve$boost_hosp, 0),
                             ve_deaths = c(ve$unvax_death, ve$prim_death, ve$boost_death, 0))
  
  #merge together
  predictions_df <- merge(covid_cases %>% group_by(age_group, boost_2_vax_status_match) %>% 
                            summarise(predicted_cases = sum(predicted_cases)),
                          ve_estimates, by.x = "boost_2_vax_status_match", by.y = "baseline_vacc", all.x = TRUE)
  
  #############
  #simulate strategies
  strat1 <- predictions_df %>% mutate(averted_cases = ifelse((boost_2_vax_status_match %in% c("Unvaccinated", "Primary Series", "Boosted (1 dose)", "Boosted (2 doses)")),
                                                             predicted_cases * ve_cases, 0) )
  
  strat1$strategy <- rep("strat1", nrow(strat1))
  strat2 <- predictions_df %>% mutate(averted_cases = ifelse((boost_2_vax_status_match %in% c("Primary Series", "Boosted (1 dose)", "Boosted (2 doses)")),
                                                             predicted_cases * ve_cases, 0))
  
  strat2$strategy <- rep("strat2", nrow(strat2))
  strat3 <- predictions_df %>% mutate(averted_cases = ifelse((boost_2_vax_status_match %in% c("Unvaccinated")),
                                                             predicted_cases * ve_cases,0))
  strat3$strategy <- rep("strat3", nrow(strat3))
  strat4 <- predictions_df %>% mutate(averted_cases = ifelse((boost_2_vax_status_match %in% c("Primary Series")),
                                                             predicted_cases * ve_cases, 0))
  strat4$strategy <- rep("strat4", nrow(strat3))
  strat5 <- predictions_df %>% mutate(averted_cases = ifelse(((age_group %in% c( "75-84 years", "85+ years")) & (boost_2_vax_status_match != "Unvaccinated")),
                                                             predicted_cases * ve_cases,0))
  strat5$strategy <- rep("strat5", nrow(strat3))
  strat6 <- predictions_df %>% mutate(averted_cases = ifelse(((age_group %in% c( "65-74 years", "75-84 years", "85+ years")) & (boost_2_vax_status_match != "Unvaccinated")),
                                                             predicted_cases * ve_cases,0))
  strat6$strategy <- rep("strat6", nrow(strat1))
  strat7 <- predictions_df %>% mutate(averted_cases = ifelse(((age_group %in% c( "50-64 years", "65-74 years", "75-84 years", "85+ years")) & (boost_2_vax_status_match != "Unvaccinated")),
                                                             predicted_cases * ve_cases, 0))
  strat7$strategy <- rep("strat7", nrow(strat1))
  
  strat8 <- predictions_df %>% mutate(averted_cases = ifelse(((age_group %in% c( "50-64 years", "65-74 years", "75-84 years", "85+ years")) & (boost_2_vax_status_match == "Unvaccinated")),
                                                             predicted_cases * ve_cases, 0))
  strat8$strategy <- rep("strat8", nrow(strat8))
  
  strat11 <- predictions_df %>% mutate(averted_cases = ifelse(((age_group %in% c( "50-64 years", "65-74 years", "75-84 years", "85+ years")) & (boost_2_vax_status_match == "Primary Series")),
                                                              predicted_cases * ve_cases, 0))
  strat11$strategy <- rep("strat11", nrow(strat11))
  
  strat14 <- predictions_df %>% mutate(averted_cases = ifelse(((age_group %in% c( "50-64 years", "65-74 years", "75-84 years", "85+ years")) & (boost_2_vax_status_match %in% c("Boosted (1 dose)", "Boosted (2 doses)"))),
                                                              predicted_cases * ve_cases, 0))
  strat14$strategy <- rep("strat14", nrow(strat14))
  
  
  strat9 <- predictions_df %>% mutate(averted_cases = ifelse(((age_group %in% c("65-74 years", "75-84 years", "85+ years")) & (boost_2_vax_status_match == "Unvaccinated")),
                                                             predicted_cases * ve_cases, 0))
  strat9$strategy <- rep("strat9", nrow(strat9))
  
  strat12 <- predictions_df %>% mutate(averted_cases = ifelse(((age_group %in% c( "65-74 years", "75-84 years", "85+ years")) & (boost_2_vax_status_match == "Primary Series")),
                                                              predicted_cases * ve_cases, 0))
  strat12$strategy <- rep("strat12", nrow(strat12))
  
  strat15 <- predictions_df %>% mutate(averted_cases = ifelse(((age_group %in% c( "65-74 years", "75-84 years", "85+ years")) & (boost_2_vax_status_match %in% c("Boosted (1 dose)", "Boosted (2 doses)"))),
                                                              predicted_cases * ve_cases, 0))
  strat15$strategy <- rep("strat15", nrow(strat15))
  
  
  strat10 <- predictions_df %>% mutate(averted_cases = ifelse(((age_group %in% c( "75-84 years", "85+ years")) & (boost_2_vax_status_match == "Unvaccinated")),
                                                              predicted_cases * ve_cases, 0))
  strat10$strategy <- rep("strat10", nrow(strat10))
  
  strat13 <- predictions_df %>% mutate(averted_cases = ifelse(((age_group %in% c( "75-84 years", "85+ years")) & (boost_2_vax_status_match == "Primary Series")),
                                                              predicted_cases * ve_cases, 0))
  strat13$strategy <- rep("strat13", nrow(strat13))
  
  strat16 <- predictions_df %>% mutate(averted_cases = ifelse(((age_group %in% c( "75-84 years", "85+ years")) & (boost_2_vax_status_match %in% c("Boosted (1 dose)", "Boosted (2 doses)"))),
                                                              predicted_cases * ve_cases, 0))
  strat16$strategy <- rep("strat16", nrow(strat16))
  
  combined_strat <- bind_rows(strat1, strat2, strat3, strat4, strat5, strat6, strat7, strat8, strat9, strat10, strat11, strat12, strat13, strat14, strat15, strat16)
  #############
  
  combined_strat_outcomes <- combined_strat %>% group_by(strategy) %>% summarise(averted_outcomes = ceiling(sum(averted_cases))) 
  
  averted_outcomes <- melt(combined_strat_outcomes %>% select(strategy, averted_outcomes)) %>% rename(averted_outcomes = value)
  
  merged_data <- merge(averted_outcomes, bivalent_booster_data %>% filter(variable == "bivalent_doses"), by = c("strategy"), all.x = TRUE) %>%
    mutate(predicted_outcomes_baseline = 1107457,
           predicted_outcomes_with_vaccine = predicted_outcomes_baseline - averted_outcomes,
           averted_outcome_by_dose = round(averted_outcomes/total_doses * 100000, 2),
           ARR = averted_outcomes/total_doses,
           NNT = 1/ARR,
           perc_averted = round(averted_outcomes/predicted_outcomes_baseline * 100, 1),
           strategy = factor(strategy, levels = c("strat1", "strat2", "strat3", "strat4", "strat5", "strat6", "strat7", "strat8", "strat9", "strat10", "strat11", "strat12", "strat13", "strat14", "strat15", "strat16"))) %>%
    arrange(strategy)
  
  vax_cases_results[i,] <- c(merged_data$NNT, merged_data$perc_averted, merged_data$averted_outcomes)
}

saveRDS(vax_cases_results, "data/results-vax-cases.RDS")
####################################################
#Paxlovid Simulations

#Create df to store results
pax_hosp_results <- data.frame(matrix(nrow = 1000, ncol = 3*11))
colnames(pax_hosp_results) <- c(sprintf("strat%d_nnt", 1:11),
                                sprintf("strat%d_perc", 1:11),
                                sprintf("strat%d_total", 1:11))
pax_death_results <- data.frame(matrix(nrow = 1000, ncol = 3*11))
colnames(pax_death_results) <- c(sprintf("strat%d_nnt", 1:11),
                                 sprintf("strat%d_perc", 1:11),
                                 sprintf("strat%d_total", 1:11))

#Run Monte Carlo simulation (N = 1000)
for (i in c(1:n)) {
  
  print(i)
  
  #create df with predictions
  covid_death <- base_df %>%  mutate(vacc_status = ifelse(boost_2_vax_status_match %in% c("Primary Series", "Boosted (1 dose)", "Boosted (2 doses)"),
                                                          "Vaccinated",
                                                          "Unvaccinated"))
  covid_death$predicted_death <- simulated_death[i,]
  covid_death$predicted_cases <- simulated_cases[i,]
  
  #create pe parameter df
  pe <- simulated_pax_params[i,]
  pe_estimates <- data.frame(baseline_vacc = c("Unvaccinated", "Vaccinated"),
                             ve_hosps = c(pe$unvax_hosp, pe$vax_hosp),
                             ve_deaths = c(pe$unvax_death, pe$vax_death))
  
  #merge together
  predictions_pax_df <- merge(covid_death 
                              %>% group_by(age_group, boost_2_vax_status_match, vacc_status) %>% 
                                summarise(predicted_cases = sum(predicted_cases),
                                          predicted_death = sum(predicted_death)),
                              pe_estimates, by.x = "vacc_status", by.y = "baseline_vacc", all.x = TRUE)
  
  
  
  #############
  #simulate strategies
  #Strategy 1: Giving Paxlovid to those 18+ years (everyone)
  strat1 <- predictions_pax_df %>% mutate(averted_death = case_when((age_group %in% c("18-49 years", "50-64 years", "65-74 years", "75-84 years", "85+ years"))  ~ predicted_death * ve_deaths,
                                                                    TRUE ~ 0),
                                          target_group_cases = case_when((age_group %in% c("18-49 years", "50-64 years", "65-74 years", "75-84 years", "85+ years")) ~ predicted_cases,
                                                                         TRUE ~ 0))
  strat1$strategy <- rep("strat1", nrow(strat1))
  
  #Strategy 2: Giving Paxlovid to those 50+ years (everyone)
  strat2 <- predictions_pax_df %>% mutate(averted_death = case_when((age_group %in% c("50-64 years", "65-74 years", "75-84 years", "85+ years"))  ~ predicted_death * ve_deaths,
                                                                    TRUE ~ 0),
                                          target_group_cases = case_when((age_group %in% c("50-64 years", "65-74 years", "75-84 years", "85+ years")) ~ predicted_cases,
                                                                         TRUE ~ 0))
  strat2$strategy <- rep("strat2", nrow(strat2))
  
  #Strategy 3: Giving Paxlovid to those 50+ years (unvaccinated and comorbidities)
  strat3 <- predictions_pax_df %>% mutate(death_unvaccinated = ifelse(age_group == "50-64 years" & vacc_status %in% c("Unvaccinated"), predicted_death, 0),
                                          death_comorbidities = ifelse(age_group == "50-64 years", predicted_death * 0.33, 0),
                                          death_both = ifelse(age_group == "50-64 years" & vacc_status %in% c("Unvaccinated"), predicted_death * 0.33, 0),
                                          death_either = death_unvaccinated + death_comorbidities - death_both,
                                          cases_unvaccinated = ifelse(age_group == "50-64 years" & vacc_status %in% c("Unvaccinated"), predicted_cases, 0),
                                          cases_comorbidities = ifelse(age_group == "50-64 years", predicted_cases * 0.33, 0),
                                          cases_both = ifelse(age_group == "50-64 years" & vacc_status %in% c("Unvaccinated"), predicted_cases * 0.33, 0),
                                          cases_either = cases_unvaccinated + cases_comorbidities - cases_both,
                                          averted_death = case_when((age_group == "50-64 years") ~ death_either * ve_deaths,
                                                                    (age_group %in% c("65-74 years", "75-84 years", "85+ years"))  ~ predicted_death * ve_deaths,
                                                                    TRUE ~ 0),
                                          target_group_cases = case_when((age_group == "50-64 years") ~ cases_either,
                                                                         (age_group %in% c("65-74 years", "75-84 years", "85+ years")) ~ predicted_cases,
                                                                         TRUE ~ 0))
  
  strat3$strategy <- rep("strat3", nrow(strat3))
  
  
  #Strategy 4: Giving Paxlovid to the those 65+ years
  strat4 <- predictions_pax_df %>% mutate(averted_death = ifelse((age_group %in% c("65-74 years", "75-84 years", "85+ years")),
                                                                 predicted_death * ve_deaths, 0),
                                          target_group_cases = ifelse((age_group %in% c("65-74 years", "75-84 years", "85+ years")),
                                                                      predicted_cases, 0))
  
  strat4$strategy <- rep("strat4", nrow(strat4))
  
  #Strategy 5: Giving Paxlovid to the those 75+ years
  strat5 <- predictions_pax_df %>% mutate(averted_death = ifelse((age_group %in% c("75-84 years", "85+ years")),
                                                                 predicted_death * ve_deaths, 0),
                                          target_group_cases = ifelse((age_group %in% c("75-84 years", "85+ years")),
                                                                      predicted_cases, 0) )
  
  strat5$strategy <- rep("strat5", nrow(strat5))  
  
  #Strategy 8: Giving Paxlovid to those Unvaccinated (65+ years)
  strat8 <- predictions_pax_df %>% mutate(averted_death = ifelse((age_group %in% c("65-74 years", "75-84 years", "85+ years") & boost_2_vax_status_match %in% c("Unvaccinated")), predicted_death * ve_deaths, 0),
                                          target_group_cases = ifelse((age_group %in% c("65-74 years", "75-84 years", "85+ years") & boost_2_vax_status_match %in% c("Unvaccinated")), predicted_cases, 0))
  
  strat8$strategy <- rep("strat8", nrow(strat8))
  
  #Strategy 10: Giving Paxlovid to those Primary Series only (65+ years)
  strat10 <- predictions_pax_df %>% mutate(averted_death = ifelse((age_group %in% c("65-74 years", "75-84 years", "85+ years") & boost_2_vax_status_match %in% c("Primary Series")), predicted_death * ve_deaths, 0),
                                           target_group_cases = ifelse((age_group %in% c("65-74 years", "75-84 years", "85+ years") & boost_2_vax_status_match %in% c("Primary Series")), predicted_cases, 0))
  
  strat10$strategy <- rep("strat10", nrow(strat10))
  
  #Strategy 12: Giving Paxlovid to those Boosted only (65+ years)
  strat12 <- predictions_pax_df %>% mutate(averted_death = ifelse((age_group %in% c("65-74 years", "75-84 years", "85+ years") & boost_2_vax_status_match %in% c("Boosted (1 dose)", "Boosted (2 doses)")), predicted_death * ve_deaths, 0),
                                           target_group_cases = ifelse((age_group %in% c("65-74 years", "75-84 years", "85+ years") & boost_2_vax_status_match %in% c("Boosted (1 dose)", "Boosted (2 doses)")), predicted_cases, 0))
  
  strat12$strategy <- rep("strat12", nrow(strat12))
  
  
  #Strategy 9: Giving Paxlovid to those Unvaccinated (75+ years)
  strat9 <- predictions_pax_df %>% mutate(averted_death = ifelse((age_group %in% c("75-84 years", "85+ years") & boost_2_vax_status_match %in% c("Unvaccinated")), predicted_death * ve_deaths, 0),
                                          target_group_cases = ifelse((age_group %in% c("75-84 years", "85+ years") & boost_2_vax_status_match %in% c("Unvaccinated")), predicted_cases, 0))
  
  strat9$strategy <- rep("strat9", nrow(strat9))
  
  #Strategy 11: Giving Paxlovid to those Primary Series only (75+ years)
  strat11 <- predictions_pax_df %>% mutate(averted_death = ifelse((age_group %in% c("75-84 years", "85+ years") & boost_2_vax_status_match %in% c("Primary Series")), predicted_death * ve_deaths, 0),
                                           target_group_cases = ifelse((age_group %in% c("75-84 years", "85+ years") & boost_2_vax_status_match %in% c("Primary Series")), predicted_cases, 0))
  
  strat11$strategy <- rep("strat11", nrow(strat11))
  
  #Strategy 13: Giving Paxlovid to those Boosted only (75+ years)
  strat13 <- predictions_pax_df %>% mutate(averted_death = ifelse((age_group %in% c("75-84 years", "85+ years") & boost_2_vax_status_match %in% c("Boosted (1 dose)", "Boosted (2 doses)")), predicted_death * ve_deaths, 0),
                                           target_group_cases = ifelse((age_group %in% c("75-84 years", "85+ years") & boost_2_vax_status_match %in% c("Boosted (1 dose)", "Boosted (2 doses)")), predicted_cases, 0))
  
  strat13$strategy <- rep("strat13", nrow(strat13))
  #############
  
  combined_strat <- bind_rows(strat1, strat2, strat3, strat4, strat5, strat8, strat9, strat10, strat11, strat12, strat13)
  
  
  combined_strat_outcomes <- combined_strat %>% group_by(strategy) %>% summarise(averted_outcomes_max = ceiling(sum(averted_death)),
                                                                                 total_death = ceiling(sum(predicted_death)),
                                                                                 target_group_cases = ceiling(sum(target_group_cases))) %>%
    mutate(averted_outcomes_cascade = ceiling(averted_outcomes_max * 0.8 * 0.895),
           treated_cascade = ceiling(target_group_cases * 0.8 *0.895),
           current_averted = ceiling(averted_outcomes_cascade * 0.3),
           current_treated = ceiling(treated_cascade * 0.3),
           actual_averted = averted_outcomes_cascade - current_averted,
           actual_treated = treated_cascade - current_treated,
           NNT = ceiling(actual_treated/actual_averted),
           perc_averted = round(actual_averted / total_death * 100, 1),
           strategy = factor(strategy, levels = c("strat1", "strat2", "strat3", "strat4", "strat5", "strat8", "strat9", "strat10", "strat11", "strat12", "strat13"))) %>%
    arrange(strategy)
  
  pax_death_results[i,] <- c(combined_strat_outcomes$NNT, combined_strat_outcomes$perc_averted, combined_strat_outcomes$actual_averted)
}

saveRDS(pax_death_results, "data/results-pax-death-new.RDS")

################################################################
#Cross-Validation Vaccine Simulations - Sensitivity Analysis

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
  covid_cases$predicted_cases <- crossval_simulated_cases[i,]
  total_cases <- ceiling(sum(covid_cases$predicted_cases))
  
  
  #merge together
  predictions_df <- covid_cases %>% group_by(age_group, boost_2_vax_status_match) %>% 
    summarise(predicted_cases = sum(predicted_cases))
  
  #############
  #simulate strategies
  strat1 <- predictions_df %>% mutate(averted_cases = ifelse((boost_2_vax_status_match %in% c("Unvaccinated", "Primary Series", "Boosted (1 dose)", "Boosted (2 doses)")),
                                                             predicted_cases, 0) )
  
  strat1$strategy <- rep("strat1", nrow(strat1))
  strat2 <- predictions_df %>% mutate(averted_cases = ifelse((boost_2_vax_status_match %in% c("Primary Series", "Boosted (1 dose)", "Boosted (2 doses)")),
                                                             predicted_cases , 0))
  
  strat2$strategy <- rep("strat2", nrow(strat2))
  strat3 <- predictions_df %>% mutate(averted_cases = ifelse((boost_2_vax_status_match %in% c("Unvaccinated")),
                                                             predicted_cases, 0))
  strat3$strategy <- rep("strat3", nrow(strat3))
  strat4 <- predictions_df %>% mutate(averted_cases = ifelse((boost_2_vax_status_match %in% c("Primary Series")),
                                                             predicted_cases, 0))
  strat4$strategy <- rep("strat4", nrow(strat3))
  strat5 <- predictions_df %>% mutate(averted_cases = ifelse(((age_group %in% c( "75-84 years", "85+ years")) & (boost_2_vax_status_match != "Unvaccinated")),
                                                             predicted_cases,0))
  strat5$strategy <- rep("strat5", nrow(strat3))
  strat6 <- predictions_df %>% mutate(averted_cases = ifelse(((age_group %in% c( "65-74 years", "75-84 years", "85+ years")) & (boost_2_vax_status_match != "Unvaccinated")),
                                                             predicted_cases,0))
  strat6$strategy <- rep("strat6", nrow(strat1))
  strat7 <- predictions_df %>% mutate(averted_cases = ifelse(((age_group %in% c( "50-64 years", "65-74 years", "75-84 years", "85+ years")) & (boost_2_vax_status_match != "Unvaccinated")),
                                                             predicted_cases, 0))
  strat7$strategy <- rep("strat7", nrow(strat1))
  
  combined_strat <- bind_rows(strat1, strat2, strat3, strat4, strat5, strat6, strat7)
  #############
  
  combined_strat_outcomes <- combined_strat %>% group_by(strategy) %>% summarise(averted_outcomes = ceiling(sum(averted_cases))) 
  
  averted_outcomes <- melt(combined_strat_outcomes %>% select(strategy, averted_outcomes)) %>% rename(averted_outcomes = value)
  
  merged_data <- merge(averted_outcomes, bivalent_booster_data %>% filter(variable == "bivalent_doses"), by = c("strategy"), all.X = TRUE) %>%
    mutate(predicted_outcomes_baseline = total_cases,
           predicted_outcomes_with_vaccine = predicted_outcomes_baseline - averted_outcomes,
           averted_outcome_by_dose = round(averted_outcomes/total_doses * 100000, 2),
           ARR = averted_outcomes/total_doses,
           NNT = 1/ARR,
           perc_averted = round(averted_outcomes/predicted_outcomes_baseline * 100, 2),
           strategy = factor(strategy, levels = c("strat1", "strat2", "strat3", "strat4", "strat5", "strat6", "strat7"))) %>%
    arrange(strategy)
  
  vax_cases_results[i,] <- c(merged_data$NNT, merged_data$perc_averted, merged_data$averted_outcomes)
}

saveRDS(vax_cases_results, "data/results-crossval-vax-cases.RDS")


#################################################################
#Waning Vaccine Simulations - Sensitivity Analysis

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
  predictions_df <- covid_cases %>% mutate(ve_cases = case_when((weeks_since_july2022 %in% c(0:4) & boost_2_vax_status_match == 'Unvaccinated') ~ ve$unvax_0_to_4,
                                                                (weeks_since_july2022 %in% c(5:26) & boost_2_vax_status_match == 'Unvaccinated') ~ ve$unvax_5_plus,
                                                                (weeks_since_july2022 %in% c(0:4) & boost_2_vax_status_match == 'Primary Series') ~ ve$prim_0_to_4,
                                                                (weeks_since_july2022 %in% c(5:9) & boost_2_vax_status_match == 'Primary Series') ~ ve$prim_5_to_9,
                                                                (weeks_since_july2022 %in% c(10:13) & boost_2_vax_status_match == 'Primary Series') ~ ve$prim_10_to_13,
                                                                (weeks_since_july2022 %in% c(14:17) & boost_2_vax_status_match == 'Primary Series') ~ ve$prim_14_to_17,
                                                                (weeks_since_july2022 %in% c(18:26) & boost_2_vax_status_match == 'Primary Series') ~ ve$prim_18_plus,
                                                                (weeks_since_july2022 %in% c(0:5) & boost_2_vax_status_match == 'Boosted (1 dose)') ~ ve$boost_0_to_5,
                                                                (weeks_since_july2022 %in% c(6:15) & boost_2_vax_status_match == 'Boosted (1 dose)') ~ ve$boost_6_to_15,
                                                                (weeks_since_july2022 %in% c(16:26) & boost_2_vax_status_match == 'Boosted (1 dose)') ~ ve$boost_16_plus,
                                                                TRUE ~ 0))
  
  #############
  #simulate strategies
  strat1 <- predictions_df %>% mutate(averted_cases = ifelse((boost_2_vax_status_match %in% c("Unvaccinated", "Primary Series", "Boosted (1 dose)", "Boosted (2 doses)")),
                                                             predicted_cases * ve_cases, 0) )
  
  strat1$strategy <- rep("strat1", nrow(strat1))
  strat2 <- predictions_df %>% mutate(averted_cases = ifelse((boost_2_vax_status_match %in% c("Primary Series", "Boosted (1 dose)", "Boosted (2 doses)")),
                                                             predicted_cases * ve_cases, 0))
  
  strat2$strategy <- rep("strat2", nrow(strat2))
  strat3 <- predictions_df %>% mutate(averted_cases = ifelse((boost_2_vax_status_match %in% c("Unvaccinated")),
                                                             predicted_cases * ve_cases,0))
  strat3$strategy <- rep("strat3", nrow(strat3))
  strat4 <- predictions_df %>% mutate(averted_cases = ifelse((boost_2_vax_status_match %in% c("Primary Series")),
                                                             predicted_cases * ve_cases, 0))
  strat4$strategy <- rep("strat4", nrow(strat3))
  strat5 <- predictions_df %>% mutate(averted_cases = ifelse(((age_group %in% c( "75-84 years", "85+ years")) & (boost_2_vax_status_match != "Unvaccinated")),
                                                             predicted_cases * ve_cases,0))
  strat5$strategy <- rep("strat5", nrow(strat3))
  strat6 <- predictions_df %>% mutate(averted_cases = ifelse(((age_group %in% c( "65-74 years", "75-84 years", "85+ years")) & (boost_2_vax_status_match != "Unvaccinated")),
                                                             predicted_cases * ve_cases,0))
  strat6$strategy <- rep("strat6", nrow(strat1))
  strat7 <- predictions_df %>% mutate(averted_cases = ifelse(((age_group %in% c( "50-64 years", "65-74 years", "75-84 years", "85+ years")) & (boost_2_vax_status_match != "Unvaccinated")),
                                                             predicted_cases * ve_cases, 0))
  strat7$strategy <- rep("strat7", nrow(strat1))
  
  combined_strat <- bind_rows(strat1, strat2, strat3, strat4, strat5, strat6, strat7)
  #############
  
  combined_strat_outcomes <- combined_strat %>% group_by(strategy) %>% summarise(averted_outcomes = ceiling(sum(averted_cases))) 
  
  averted_outcomes <- melt(combined_strat_outcomes %>% select(strategy, averted_outcomes)) %>% rename(averted_outcomes = value)
  
  merged_data <- merge(averted_outcomes, bivalent_booster_data %>% filter(variable == "bivalent_doses"), by = c("strategy"), all.X = TRUE) %>%
    mutate(predicted_outcomes_baseline = 1107457,
           predicted_outcomes_with_vaccine = predicted_outcomes_baseline - averted_outcomes,
           averted_outcome_by_dose = round(averted_outcomes/total_doses * 100000, 2),
           ARR = averted_outcomes/total_doses,
           NNT = 1/ARR,
           perc_averted = round(averted_outcomes/predicted_outcomes_baseline * 100, 2),
           strategy = factor(strategy, levels = c("strat1", "strat2", "strat3", "strat4", "strat5", "strat6", "strat7"))) %>%
    arrange(strategy)
  
  vax_cases_results[i,] <- c(merged_data$NNT, merged_data$perc_averted, merged_data$averted_outcomes)
}

saveRDS(vax_cases_results, "data/results-waning-vax-cases.RDS")

#################################################################
#Prospective Vaccination Simulations - Sensitivity Analysis

#Create df to store results
vax_cases_results <- data.frame(matrix(nrow = 1000, ncol = 3*7))
colnames(vax_cases_results) <- c(sprintf("strat%d_nnt", 1:7),
                                 sprintf("strat%d_perc", 1:7),
                                 sprintf("strat%d_total", 1:7))
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
  covid_cases <- base_df 
  covid_cases$predicted_cases <- simulated_cases[i,]
  
  #create ve parameter df
  ve <- simulated_vax_params[i,]
  ve_estimates <- data.frame(baseline_vacc = c("Unvaccinated", "Primary Series", "Boosted (1 dose)","Boosted (2 doses)"),
                             ve_cases = c(ve$unvax_case, ve$prim_case, ve$boost_case, 0),
                             ve_hosps = c(ve$unvax_hosp, ve$prim_hosp, ve$boost_hosp, 0),
                             ve_deaths = c(ve$unvax_death, ve$prim_death, ve$boost_death, 0))
  
  #merge together
  predictions_df <- merge(covid_cases %>% 
                            mutate(vax_coverage = case_when((weeks_since_july2022 %in% c(0:4)) ~ 0,
                                                            (weeks_since_july2022 %in% c(5:8)) ~ 0.2,
                                                            (weeks_since_july2022 %in% c(9:13)) ~ 0.4,
                                                            (weeks_since_july2022 %in% c(14:18)) ~ 0.6,
                                                            (weeks_since_july2022 %in% c(19:22)) ~ 0.8,
                                                            (weeks_since_july2022 %in% c(23:26)) ~ 1,
                                                            TRUE ~ 0)),
                          ve_estimates, by.x = "boost_2_vax_status_match", by.y = "baseline_vacc", all.x = TRUE)
  
  #############
  #simulate strategies
  strat1 <- predictions_df %>% mutate(averted_cases = ifelse((boost_2_vax_status_match %in% c("Unvaccinated", "Primary Series", "Boosted (1 dose)", "Boosted (2 doses)")),
                                                             predicted_cases * ve_cases * vax_coverage, 0) )
  
  strat1$strategy <- rep("strat1", nrow(strat1))
  strat2 <- predictions_df %>% mutate(averted_cases = ifelse((boost_2_vax_status_match %in% c("Primary Series", "Boosted (1 dose)", "Boosted (2 doses)")),
                                                             predicted_cases * ve_cases * vax_coverage, 0))
  
  strat2$strategy <- rep("strat2", nrow(strat2))
  strat3 <- predictions_df %>% mutate(averted_cases = ifelse((boost_2_vax_status_match %in% c("Unvaccinated")),
                                                             predicted_cases * ve_cases * vax_coverage,0))
  strat3$strategy <- rep("strat3", nrow(strat3))
  strat4 <- predictions_df %>% mutate(averted_cases = ifelse((boost_2_vax_status_match %in% c("Primary Series")),
                                                             predicted_cases * ve_cases * vax_coverage, 0))
  strat4$strategy <- rep("strat4", nrow(strat3))
  strat5 <- predictions_df %>% mutate(averted_cases = ifelse(((age_group %in% c( "75-84 years", "85+ years")) & (boost_2_vax_status_match != "Unvaccinated")),
                                                             predicted_cases * ve_cases * vax_coverage,0))
  strat5$strategy <- rep("strat5", nrow(strat3))
  strat6 <- predictions_df %>% mutate(averted_cases = ifelse(((age_group %in% c( "65-74 years", "75-84 years", "85+ years")) & (boost_2_vax_status_match != "Unvaccinated")),
                                                             predicted_cases * ve_cases * vax_coverage,0))
  strat6$strategy <- rep("strat6", nrow(strat1))
  strat7 <- predictions_df %>% mutate(averted_cases = ifelse(((age_group %in% c( "50-64 years", "65-74 years", "75-84 years", "85+ years")) & (boost_2_vax_status_match != "Unvaccinated")),
                                                             predicted_cases * ve_cases * vax_coverage, 0))
  strat7$strategy <- rep("strat7", nrow(strat1))
  
  combined_strat <- bind_rows(strat1, strat2, strat3, strat4, strat5, strat6, strat7)
  
  combined_strat_outcomes <- combined_strat %>% group_by(strategy) %>% summarise(averted_outcomes = ceiling(sum(averted_cases))) 
  
  averted_outcomes <- melt(combined_strat_outcomes %>% select(strategy, averted_outcomes)) %>% rename(averted_outcomes = value)
  
  merged_data <- merge(averted_outcomes, bivalent_booster_data %>% filter(variable == "bivalent_doses"), by = c("strategy"), all.X = TRUE) %>%
    mutate(predicted_outcomes_baseline = 1107457,
           predicted_outcomes_with_vaccine = predicted_outcomes_baseline - averted_outcomes,
           averted_outcome_by_dose = round(averted_outcomes/total_doses * 100000, 2),
           ARR = averted_outcomes/total_doses,
           NNT = ceiling(1/ARR),
           perc_averted = round(averted_outcomes/predicted_outcomes_baseline * 100, 1),
           strategy = factor(strategy, levels = c("strat1", "strat2", "strat3", "strat4", "strat5", "strat6", "strat7"))) %>%
    arrange(strategy)
  
  vax_cases_results[i,] <- c(merged_data$NNT, merged_data$perc_averted, merged_data$averted_outcomes)
}

saveRDS(vax_cases_results, "data/results-prospective-vax-death.RDS")

#################################################################
#Evaluate the 95% UIs for the results
results_pax_death <- readRDS("data/results-pax-death-new.RDS")
results_pax_hosp <- readRDS("data/results-pax-hosp-new.RDS")
results_vax_hosp <- readRDS("data/results-vax-hosp.RDS")
results_vax_cases <- readRDS("data/results-vax-cases.RDS")
results_vax_death <- readRDS("data/results-vax-death.RDS")
results_waning_vax_cases <- readRDS("data/results-waning-vax-cases.RDS")
results_prospective_vax_hosp <- readRDS("data/results-prospective-vax-hosp.RDS")
results_crossval_vax_cases <- readRDS("data/results-crossval-vax-cases.RDS")
results_prospective_vax_death <- readRDS("data/results-prospective-vax-death.RDS")



inspection <- results_pax_death %>% apply(2, quantile, probs=c(0.025, 0.975))






