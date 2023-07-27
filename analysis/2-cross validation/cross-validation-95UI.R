###################################################################################################
#Title: Cross-Validation (with 95% UI)
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
base_df_crossval <- readRDS("data/base-dataframe-crossval-hosp.RDS") 
crossval_simulated_cases <- readRDS('data/simulated-cases-crossval.RDS')
crossval_simulated_hosp <- readRDS('data/simulated-hosp-crossval.RDS')
crossval_simulated_death <- readRDS('data/simulated-death-crossval.RDS')

#Number of simulations
n <- 1000



#Cross-Validation Simulations

#Create df to store results
vax_cases_results <- data.frame(matrix(nrow = 1000, ncol = 1*7)) 
colnames(vax_cases_results) <- c(sprintf("strat%d_predicted_outcomes", 1:7))
vax_hosp_results <- data.frame(matrix(nrow = 1000, ncol = 1*7)) 
colnames(vax_hosp_results) <- c(sprintf("strat%d_predicted_outcomes", 1:7))
vax_death_results <- data.frame(matrix(nrow = 1000, ncol = 1*7)) 
colnames(vax_death_results) <- c(sprintf("strat%d_predicted_outcomes", 1:7))

#Run Monte Carlo simulation (N = 1000)
#NOTE: Switch out df names for respective outcome analyses.
for (i in c(1:n)) {
  print(i)
  #create df with predictions
  covid_df <- base_df_crossval
  covid_df$predicted_outcomes <- crossval_simulated_death[i,] #CHANGE HERE ('crossval_simulated_cases', 'crossval_simulated_hosp', 'crossval_simulated_death')
  total_outcomes <- ceiling(sum(covid_df$predicted_outcomes))
  
  
  #merge together
  predictions_df <- covid_df %>% group_by(age_group, boost_2_vax_status_match) %>% 
    summarise(predicted_outcomes = sum(predicted_outcomes))
  
  #############
  #simulate strategies
  strat1 <- predictions_df %>% mutate(predicted_outcomes = ifelse((boost_2_vax_status_match %in% c("Unvaccinated", "Primary Series", "Boosted (1 dose)", "Boosted (2 doses)")),
                                                                  predicted_outcomes, 0) )
  
  strat1$strategy <- rep("strat1", nrow(strat1))
  strat2 <- predictions_df %>% mutate(predicted_outcomes = ifelse((boost_2_vax_status_match %in% c("Primary Series", "Boosted (1 dose)", "Boosted (2 doses)")),
                                                                  predicted_outcomes , 0))
  
  strat2$strategy <- rep("strat2", nrow(strat2))
  strat3 <- predictions_df %>% mutate(predicted_outcomes = ifelse((boost_2_vax_status_match %in% c("Unvaccinated")),
                                                                  predicted_outcomes, 0))
  strat3$strategy <- rep("strat3", nrow(strat3))
  strat4 <- predictions_df %>% mutate(predicted_outcomes = ifelse((boost_2_vax_status_match %in% c("Primary Series")),
                                                                  predicted_outcomes, 0))
  strat4$strategy <- rep("strat4", nrow(strat3))
  strat5 <- predictions_df %>% mutate(predicted_outcomes = ifelse(((age_group %in% c( "75-84 years", "85+ years")) & (boost_2_vax_status_match != "Unvaccinated")),
                                                                  predicted_outcomes,0))
  strat5$strategy <- rep("strat5", nrow(strat3))
  strat6 <- predictions_df %>% mutate(predicted_outcomes = ifelse(((age_group %in% c( "65-74 years", "75-84 years", "85+ years")) & (boost_2_vax_status_match != "Unvaccinated")),
                                                                  predicted_outcomes,0))
  strat6$strategy <- rep("strat6", nrow(strat1))
  strat7 <- predictions_df %>% mutate(predicted_outcomes = ifelse(((age_group %in% c( "50-64 years", "65-74 years", "75-84 years", "85+ years")) & (boost_2_vax_status_match != "Unvaccinated")),
                                                                  predicted_outcomes, 0))
  strat7$strategy <- rep("strat7", nrow(strat1))
  
  combined_strat <- bind_rows(strat1, strat2, strat3, strat4, strat5, strat6, strat7)
  #############
  
  combined_strat_outcomes <- combined_strat %>% group_by(strategy) %>% summarise(predicted_outcomes = ceiling(sum(predicted_outcomes))) 
  
  vax_death_results[i,] <- c(combined_strat_outcomes$predicted_outcomes) #CHANGE HERE ('vax_death_results', 'vax_hosp_results', 'vax_cases_results')
}

saveRDS(vax_death_results, "data/results-crossval-death.RDS")  #CHANGE HERE ('vax_death_results', 'vax_hosp_results', 'vax_cases_results')




results_crossval_cases <- readRDS("data/results-crossval-cases.RDS")
results_crossval_hosp <- readRDS("data/results-crossval-hosp.RDS")
results_crossval_death <- readRDS("data/results-crossval-death.RDS")


inspection <- results_crossval_death %>% apply(2, quantile, probs=c(0.025, 0.975))
