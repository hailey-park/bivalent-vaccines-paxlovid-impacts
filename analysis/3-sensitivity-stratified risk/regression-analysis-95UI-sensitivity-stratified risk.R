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
bivalent_booster_data <- read.csv("outcome data clean/bivalent_booster_doses_clean_70maxuptake_v2.csv")[,-1]
simulated_vax_waning_params <- readRDS('data/simulated-ve-waning-parameters.RDS')
bivalent_coverage <- read.csv("outcome data clean/bivalent_coverage_age_vax.csv")[,-1]


#Number of simulations
n <- 1000

#####################################################################
#Bivalent Vaccine Simulations (Outcomes: Hospitalizations, Deaths) 
#Create df to store results
vax_hosp_results <- data.frame(matrix(nrow = 1000, ncol = 3*9))
colnames(vax_hosp_results) <- c(sprintf("strat%d_nnt", 8:16),
                                 sprintf("strat%d_perc", 8:16),
                                 sprintf("strat%d_total", 8:16))
vax_death_results <- data.frame(matrix(nrow = 1000, ncol = 3*9))
colnames(vax_death_results) <- c(sprintf("strat%d_nnt", 8:16),
                                sprintf("strat%d_perc", 8:16),
                                sprintf("strat%d_total", 8:16))

#Run Monte Carlo simulation (N = 1000)
#NOTE: Switch out df names for respective outcome analyses.
for (i in c(1:n)) {
  print(i)
  #create df with predictions
  covid_df <- base_df 
  covid_df$predicted_outcomes <- simulated_hosp[i,] #CHOOSE 'simulated_hosp' or 'simulated_death'
  
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
    mutate(predicted_outcomes = round(predicted_outcomes))
  
  #simulate strategies
  #make sure to switch out vaccine effectiveness variables with respective outcomes ('ve_hosps', 've_deaths')
  strat8 <- predictions_df %>% mutate(averted_diff = ifelse(((age_group %in% c( "50-64 years", "65-74 years", "75-84 years", "85+ years")) & (boost_2_vax_status_match == "Unvaccinated")),
                                                            predicted_outcomes * ve_hosps * (0.7 - bivalent_coverage), 0),
                                      total_predicted = predicted_outcomes * (1 - (bivalent_coverage * ve_hosps)))
  strat8$strategy <- rep("strat8", nrow(strat8))

  strat9 <- predictions_df %>% mutate(averted_diff = ifelse(((age_group %in% c("65-74 years", "75-84 years", "85+ years")) & (boost_2_vax_status_match == "Unvaccinated")),
                                                            predicted_outcomes * ve_hosps * (0.7 - bivalent_coverage), 0),
                                      total_predicted = predicted_outcomes * (1 - (bivalent_coverage * ve_hosps)))
  strat9$strategy <- rep("strat9", nrow(strat9))

  strat10 <- predictions_df %>% mutate(averted_diff = ifelse(((age_group %in% c( "75-84 years", "85+ years")) & (boost_2_vax_status_match == "Unvaccinated")),
                                                             predicted_outcomes * ve_hosps * (0.7 - bivalent_coverage), 0),
                                       total_predicted = predicted_outcomes * (1 - (bivalent_coverage * ve_hosps)))
  strat10$strategy <- rep("strat10", nrow(strat10))
  
  strat11 <- predictions_df %>% mutate(averted_diff = ifelse(((age_group %in% c( "50-64 years", "65-74 years", "75-84 years", "85+ years")) & (boost_2_vax_status_match == "Primary Series")),
                                                             predicted_outcomes * ve_hosps * (0.7 - bivalent_coverage), 0),
                                       total_predicted = predicted_outcomes * (1 - (bivalent_coverage * ve_hosps)))
  strat11$strategy <- rep("strat11", nrow(strat11))

  strat12 <- predictions_df %>% mutate(averted_diff = ifelse(((age_group %in% c( "65-74 years", "75-84 years", "85+ years")) & (boost_2_vax_status_match == "Primary Series")),
                                                             predicted_outcomes * ve_hosps * (0.7 - bivalent_coverage), 0),
                                       total_predicted = predicted_outcomes * (1 - (bivalent_coverage * ve_hosps)))
  strat12$strategy <- rep("strat12", nrow(strat12))

  strat13 <- predictions_df %>% mutate(averted_diff = ifelse(((age_group %in% c( "75-84 years", "85+ years")) & (boost_2_vax_status_match == "Primary Series")),
                                                             predicted_outcomes * ve_hosps * (0.7 - bivalent_coverage), 0),
                                       total_predicted = predicted_outcomes * (1 - (bivalent_coverage * ve_hosps)))
  strat13$strategy <- rep("strat13", nrow(strat13))
  
  strat14 <- predictions_df %>% mutate(averted_diff = ifelse(((age_group %in% c( "50-64 years", "65-74 years", "75-84 years", "85+ years")) & (boost_2_vax_status_match %in% c("Boosted (1 dose)", "Boosted (2 doses)"))),
                                                             predicted_outcomes * ve_hosps * (0.7 - bivalent_coverage), 0),
                                       total_predicted = predicted_outcomes * (1 - (bivalent_coverage * ve_hosps)))
  strat14$strategy <- rep("strat14", nrow(strat14))
  
  strat15 <- predictions_df %>% mutate(averted_diff = ifelse(((age_group %in% c( "65-74 years", "75-84 years", "85+ years")) & (boost_2_vax_status_match %in% c("Boosted (1 dose)", "Boosted (2 doses)"))),
                                                             predicted_outcomes * ve_hosps * (0.7 - bivalent_coverage), 0),
                                       total_predicted = predicted_outcomes * (1 - (bivalent_coverage * ve_hosps)))
  strat15$strategy <- rep("strat15", nrow(strat15))
  
  strat16 <- predictions_df %>% mutate(averted_diff = ifelse(((age_group %in% c( "75-84 years", "85+ years")) & (boost_2_vax_status_match %in% c("Boosted (1 dose)", "Boosted (2 doses)"))),
                                                             predicted_outcomes * ve_hosps * (0.7 - bivalent_coverage), 0),
                                       total_predicted = predicted_outcomes * (1 - (bivalent_coverage * ve_hosps)))
  strat16$strategy <- rep("strat16", nrow(strat16))
  
  combined_strat <- bind_rows(strat8, strat9, strat10, strat11, strat12, strat13, strat14, strat15, strat16)
  #############
  
  combined_strat_outcomes <- combined_strat %>% group_by(strategy) %>% summarise(averted_outcomes_diff = round(sum(averted_diff)),
                                                                                 total_predicted_outcomes = round(sum(total_predicted))) 
  
  merged_data <- merge(combined_strat_outcomes, bivalent_booster_data, by = c("strategy"), all.X = TRUE) %>%
    mutate(ARR = averted_outcomes_diff/(total_doses),
           NNT = ceiling(1/ARR),
           perc_averted = round(averted_outcomes_diff/total_predicted_outcomes * 100, 1),
           strategy = factor(strategy, levels = c("strat8", "strat9", "strat10", "strat11", "strat12", "strat13", "strat14", "strat15", "strat16")))%>%  
    select(-c(variable))  %>% arrange(strategy)
  
  
  vax_hosp_results[i,] <- c(merged_data$NNT, merged_data$perc_averted, merged_data$averted_outcomes_diff)
}

saveRDS(vax_hosp_results, "data/results-sensitivity-stratified risk-vax-hosp.RDS") #CHANGE HERE ('hosp' vs. 'death')
#####################################################################
#Bivalent Vaccine Simulations with Waning (Outcomes: Cases) 

#Create df to store results
vvax_cases_results <- data.frame(matrix(nrow = 1000, ncol = 3*9))
colnames(vax_cases_results) <- c(sprintf("strat%d_nnt", 8:16),
                                 sprintf("strat%d_perc", 8:16),
                                 sprintf("strat%d_total", 8:16))

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
    mutate(predicted_outcomes = round(predicted_cases))
  
  #############
  #simulate strategies
  
  strat8 <- predictions_df %>% mutate(averted_diff = ifelse(((age_group %in% c( "50-64 years", "65-74 years", "75-84 years", "85+ years")) & (boost_2_vax_status_match == "Unvaccinated")),
                                                            predicted_outcomes * ve_cases * (0.7 - bivalent_coverage), 0),
                                      total_predicted = predicted_outcomes * (1 - (bivalent_coverage * ve_cases)))
  strat8$strategy <- rep("strat8", nrow(strat8))

  strat9 <- predictions_df %>% mutate(averted_diff = ifelse(((age_group %in% c("65-74 years", "75-84 years", "85+ years")) & (boost_2_vax_status_match == "Unvaccinated")),
                                                            predicted_outcomes * ve_cases * (0.7 - bivalent_coverage), 0),
                                      total_predicted = predicted_outcomes * (1 - (bivalent_coverage * ve_cases)))
  strat9$strategy <- rep("strat9", nrow(strat9))

  strat10 <- predictions_df %>% mutate(averted_diff = ifelse(((age_group %in% c( "75-84 years", "85+ years")) & (boost_2_vax_status_match == "Unvaccinated")),
                                                             predicted_outcomes * ve_cases * (0.7 - bivalent_coverage), 0),
                                       total_predicted = predicted_outcomes * (1 - (bivalent_coverage * ve_cases)))
  strat10$strategy <- rep("strat10", nrow(strat10))
  
  strat11 <- predictions_df %>% mutate(averted_diff = ifelse(((age_group %in% c( "50-64 years", "65-74 years", "75-84 years", "85+ years")) & (boost_2_vax_status_match == "Primary Series")),
                                                             predicted_outcomes * ve_cases * (0.7 - bivalent_coverage), 0),
                                       total_predicted = predicted_outcomes * (1 - (bivalent_coverage * ve_cases)))
  strat11$strategy <- rep("strat11", nrow(strat11))

  strat12 <- predictions_df %>% mutate(averted_diff = ifelse(((age_group %in% c( "65-74 years", "75-84 years", "85+ years")) & (boost_2_vax_status_match == "Primary Series")),
                                                             predicted_outcomes * ve_cases * (0.7 - bivalent_coverage), 0),
                                       total_predicted = predicted_outcomes * (1 - (bivalent_coverage * ve_cases)))
  strat12$strategy <- rep("strat12", nrow(strat12))

  strat13 <- predictions_df %>% mutate(averted_diff = ifelse(((age_group %in% c( "75-84 years", "85+ years")) & (boost_2_vax_status_match == "Primary Series")),
                                                             predicted_outcomes * ve_cases * (0.7 - bivalent_coverage), 0),
                                       total_predicted = predicted_outcomes * (1 - (bivalent_coverage * ve_cases)))
  strat13$strategy <- rep("strat13", nrow(strat13))
  
  
  strat14 <- predictions_df %>% mutate(averted_diff = ifelse(((age_group %in% c( "50-64 years", "65-74 years", "75-84 years", "85+ years")) & (boost_2_vax_status_match %in% c("Boosted (1 dose)", "Boosted (2 doses)"))),
                                                             predicted_outcomes * ve_cases * (0.7 - bivalent_coverage), 0),
                                       total_predicted = predicted_outcomes * (1 - (bivalent_coverage * ve_cases)))
  strat14$strategy <- rep("strat14", nrow(strat14))
  
  strat15 <- predictions_df %>% mutate(averted_diff = ifelse(((age_group %in% c( "65-74 years", "75-84 years", "85+ years")) & (boost_2_vax_status_match %in% c("Boosted (1 dose)", "Boosted (2 doses)"))),
                                                             predicted_outcomes * ve_cases * (0.7 - bivalent_coverage), 0),
                                       total_predicted = predicted_outcomes * (1 - (bivalent_coverage * ve_cases)))
  strat15$strategy <- rep("strat15", nrow(strat15))
  
  strat16 <- predictions_df %>% mutate(averted_diff = ifelse(((age_group %in% c( "75-84 years", "85+ years")) & (boost_2_vax_status_match %in% c("Boosted (1 dose)", "Boosted (2 doses)"))),
                                                             predicted_outcomes * ve_cases * (0.7 - bivalent_coverage), 0),
                                       total_predicted = predicted_outcomes * (1 - (bivalent_coverage * ve_cases)))
  strat16$strategy <- rep("strat16", nrow(strat16))
  
  combined_strat <- bind_rows(strat8, strat9, strat10, strat11, strat12, strat13, strat14, strat15, strat16)
  #############

  combined_strat_outcomes <- combined_strat %>% group_by(strategy) %>% summarise(averted_outcomes_diff = round(sum(averted_diff)),
                                                                                 total_predicted_outcomes = round(sum(total_predicted))) 
  
  merged_data <- merge(combined_strat_outcomes, bivalent_booster_data, by = c("strategy"), all.X = TRUE) %>%
    mutate(ARR = averted_outcomes_diff/(total_doses),
           NNT = ceiling(1/ARR),
           perc_averted = round(averted_outcomes_diff/total_predicted_outcomes * 100, 1),
           strategy = factor(strategy, levels = c("strat8", "strat9", "strat10", "strat11", "strat12", "strat13", "strat14", "strat15", "strat16"))) %>%  
    select(-c(variable))  %>% arrange(strategy)
  
  
  vax_cases_results[i,] <- c(merged_data$NNT, merged_data$perc_averted, merged_data$averted_outcomes_diff)
}

saveRDS(vax_cases_results, "data/results-sensitivity-stratified risk-vax-cases.RDS")

#####################################################################
#Paxlovid Simulations

#Create df to store results
pax_hosp_results <- data.frame(matrix(nrow = 1000, ncol = 3*6))
colnames(pax_hosp_results) <- c(sprintf("strat%d_nnt", 6:12,
                                 sprintf("strat%d_perc", 6:12,
                                 sprintf("strat%d_total", 6:12)
pax_hosp_results <- data.frame(matrix(nrow = 1000, ncol = 3*6))
colnames(pax_hosp_results) <- c(sprintf("strat%d_nnt", 6:12,
                                sprintf("strat%d_perc", 6:12,
                                sprintf("strat%d_total", 6:12)

#Run Monte Carlo simulation (N = 1000)
#NOTE: Switch out df names for respective outcome analyses.                                   
                                        
for (i in c(1:n)) {
  
  print(i)

  #create df with predictions
  #CHANGE ('covid_hosp' or 'covid_death')
  covid_hosp <- base_df %>%  mutate(vacc_status = ifelse(boost_2_vax_status_match %in% c("Primary Series", "Boosted (1 dose)", "Boosted (2 doses)"),
                                                                                                      "Vaccinated",
                                                                                                      "Unvaccinated"))
  covid_hosp$predicted_hosp <- simulated_hosp[i,] #CHANGE ('simulated_hosp' or 'simulated_death')
  covid_hosp$predicted_cases <- simulated_cases[i,]
  
  #create pe parameter df
  pe <- simulated_pax_params[i,]
  pe_estimates <- data.frame(baseline_vacc = c("Unvaccinated", "Vaccinated"),
                             pe_hosps = c(pe$unvax_hosp, pe$vax_hosp),
                             pe_deaths = c(pe$unvax_death, pe$vax_death))
  
  #merge together and add age-specific baseline Paxlovid uptake and maximum Paxlovid coverage
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
  
  #Additional Paxlovid Strategies
  #make sure to switch out variables with respective outcomes ('hosp' vs. 'death')
  
  #Strategy 7: Giving Paxlovid to those Unvaccinated (65+ years)
  strat7 <- predictions_pax_df %>% mutate(averted_hosp_diff = ifelse((age_group %in% c("65-74 years", "75-84 years", "85+ years") & boost_2_vax_status_match %in% c("Unvaccinated")), predicted_hosp * pe_hosps * (max_pax_uptake - baseline_pax_uptake), 0),
                                          target_group_cases_diff = ifelse((age_group %in% c("65-74 years", "75-84 years", "85+ years") & boost_2_vax_status_match %in% c("Unvaccinated")), predicted_cases * (max_pax_uptake - baseline_pax_uptake), 0),
                                          total_predicted = predicted_hosp * (1 - (baseline_pax_uptake * pe_hosps)))
  
  strat7$strategy <- rep("strat7", nrow(strat7))
  
  #Strategy 8: Giving Paxlovid to those Unvaccinated (75+ years)
  strat8 <- predictions_pax_df %>% mutate(averted_hosp_diff = ifelse((age_group %in% c("75-84 years", "85+ years") & boost_2_vax_status_match %in% c("Unvaccinated")), predicted_hosp * pe_hosps * (max_pax_uptake - baseline_pax_uptake), 0),
                                          target_group_cases_diff = ifelse((age_group %in% c("75-84 years", "85+ years") & boost_2_vax_status_match %in% c("Unvaccinated")), predicted_cases * (max_pax_uptake - baseline_pax_uptake), 0),
                                          total_predicted = predicted_hosp * (1 - (baseline_pax_uptake * pe_hosps)))
  
  strat8$strategy <- rep("strat8", nrow(strat8))
  
  
  #Strategy 9:Giving Paxlovid to those Primary Series only (65+ years)
  strat9 <- predictions_pax_df %>% mutate(averted_hosp_diff = ifelse((age_group %in% c("65-74 years", "75-84 years", "85+ years") & boost_2_vax_status_match %in% c("Primary Series")), predicted_hosp * pe_hosps * (max_pax_uptake - baseline_pax_uptake), 0),
                                          target_group_cases_diff = ifelse((age_group %in% c("65-74 years", "75-84 years", "85+ years") & boost_2_vax_status_match %in% c("Primary Series")), predicted_cases * (max_pax_uptake - baseline_pax_uptake), 0),
                                          total_predicted = predicted_hosp * (1 - (baseline_pax_uptake * pe_hosps)))
  
  strat9$strategy <- rep("strat9", nrow(strat9))
  
  #Strategy 10: Giving Paxlovid to those Primary Series only (75+ years)
  strat10 <- predictions_pax_df %>% mutate(averted_hosp_diff = ifelse((age_group %in% c("75-84 years", "85+ years") & boost_2_vax_status_match %in% c("Primary Series")), predicted_hosp * pe_hosps * (max_pax_uptake - baseline_pax_uptake), 0),
                                          target_group_cases_diff = ifelse((age_group %in% c("75-84 years", "85+ years") & boost_2_vax_status_match %in% c("Primary Series")), predicted_cases * (max_pax_uptake - baseline_pax_uptake), 0),
                                          total_predicted = predicted_hosp * (1 - (baseline_pax_uptake * pe_hosps)))
  
  strat10$strategy <- rep("strat10", nrow(strat10))
  
  #Strategy 11: Giving Paxlovid to those Boosted only (65+ years)
  strat11 <- predictions_pax_df %>% mutate(averted_hosp_diff = ifelse((age_group %in% c("65-74 years", "75-84 years", "85+ years") & boost_2_vax_status_match %in% c("Boosted (1 dose)", "Boosted (2 doses)")), predicted_hosp * pe_hosps * (max_pax_uptake - baseline_pax_uptake), 0),
                                          target_group_cases_diff = ifelse((age_group %in% c("65-74 years", "75-84 years", "85+ years") & boost_2_vax_status_match %in% c("Boosted (1 dose)", "Boosted (2 doses)")), predicted_cases * (max_pax_uptake - baseline_pax_uptake), 0),
                                          total_predicted = predicted_hosp * (1 - (baseline_pax_uptake * pe_hosps)))
  
  strat11$strategy <- rep("strat11", nrow(strat11))
  
  #Strategy 12: Giving Paxlovid to those Boosted only (75+ years)
  strat12 <- predictions_pax_df %>% mutate(averted_hosp_diff = ifelse((age_group %in% c("75-84 years", "85+ years") & boost_2_vax_status_match %in% c("Boosted (1 dose)", "Boosted (2 doses)")), predicted_hosp * pe_hosps * (max_pax_uptake - baseline_pax_uptake), 0),
                                          target_group_cases_diff = ifelse((age_group %in% c("75-84 years", "85+ years") & boost_2_vax_status_match %in% c("Boosted (1 dose)", "Boosted (2 doses)")), predicted_cases * (max_pax_uptake - baseline_pax_uptake), 0),
                                          total_predicted = predicted_hosp * (1 - (baseline_pax_uptake * pe_hosps)))
  
  strat12$strategy <- rep("strat12", nrow(strat12))
  combined_strat <- bind_rows(strat7, strat8, strat9, strat10, strat11, strat12)
  
  
  #######################################
  combined_strat_outcomes <- combined_strat %>% group_by(strategy) %>% summarise(averted_outcomes_diff = ceiling(sum(averted_hosp_diff)),
                                                                                 total_hosp = ceiling(sum(total_predicted)),
                                                                                 target_group_cases_diff = ceiling(sum(target_group_cases_diff))) %>%
    mutate(averted_outcomes_cascade = ceiling(averted_outcomes_diff * 0.8 * 0.895),
           treated_cascade = ceiling(target_group_cases_diff * 0.8 * 0.895),
           NNT = ceiling(treated_cascade/averted_outcomes_cascade),
           perc_averted = round(averted_outcomes_cascade/ total_hosp * 100, 1),
           strategy = factor(strategy, levels = c("strat7", "strat8", "strat9", "strat10", "strat11", "strat12"))) %>%
    arrange(strategy)
  
  pax_hosp_results[i,] <- c(combined_strat_outcomes$NNT, combined_strat_outcomes$perc_averted, combined_strat_outcomes$averted_outcomes_cascade)
}

saveRDS(pax_hosp_results, "data/results-sensitivity-stratified risk-pax-hosp.RDS") #CHANGE ('hosp' or 'death')


#################################################################
#Evaluate the 95% UIs for the results
results_pax_death <- readRDS("data/results-sensitivity-stratified risk-pax-death.RDS")
results_pax_hosp <- readRDS("data/results-sensitivity-stratified risk-pax-hosp.RDS")
results_vax_hosp <- readRDS("data/results-sensitivity-stratified risk-vax-hosp.RDS")
results_vax_cases <- readRDS("data/results-sensitivity-stratified risk-vax-cases.RDS")
results_vax_death <- readRDS("data/results-sensitivity-stratified risk-vax-death.RDS")

inspection <- results_vax_cases %>% apply(2, quantile, probs=c(0.025, 0.975))





