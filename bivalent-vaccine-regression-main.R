###################################################################################################
#Title: Bivalent Booster Regression Analysis
#Author: Hailey Park
#Date: January 30, 2023
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
covid_case_data <- read.csv("outcome data clean/covid_data_clean_cases_crossval-6mo.csv")[,-1]
covid_hosp_data <- read.csv("outcome data clean/covid_data_clean_hosp_crossval-6mo.csv")[,-1]
covid_death_data <- read.csv("outcome data clean/covid_data_clean_death_crossval-6mo.csv")[,-1]
bivalent_booster_data <- read.csv("outcome data clean/bivalent_booster_doses_clean_new.csv")[,-1]

#Make dataframes for each outcome (cases, hosp, deaths)
covid_cases <- covid_case_data %>%
  mutate(weeks_since_july2022 = as.integer(floor(difftime(strptime(dtepisode, format="%Y-%m-%d"), 
                                                          strptime("2022-01-23", format="%Y-%m-%d"), 
                                                          units = 'weeks')))) %>%
  group_by(weeks_since_july2022, age_group, boost_2_vax_status_match, .drop = FALSE) %>%
  summarise(num_cases = n()) 

covid_hosp <- covid_hosp_data  %>%
  mutate(weeks_since_july2022 = as.integer(floor(difftime(strptime(dt_hosp, format="%Y-%m-%d"), 
                                                          strptime("2022-01-23", format="%Y-%m-%d"), 
                                                          units = 'weeks')))) %>%
  group_by(weeks_since_july2022, age_group, boost_2_vax_status_match, .drop = FALSE) %>%
  summarise(num_hosp = n()) %>% mutate(num_hosp = (num_hosp * (1/0.75)))

covid_death <- covid_death_data %>%
  mutate(weeks_since_july2022 = as.integer(floor(difftime(strptime(death_date_both_or, format="%Y-%m-%d"), 
                                                          strptime("2022-01-23", format="%Y-%m-%d"), 
                                                          units = 'weeks')))) %>%
  group_by(weeks_since_july2022, age_group, boost_2_vax_status_match, .drop = FALSE) %>%
  summarise(num_death = n()) 

#Create dataframes for VE estimates for averting outcomes
ve_estimates <- data.frame(baseline_vacc = c("Unvaccinated", "Primary Series", "Boosted (1 dose)","Boosted (2 doses)"),
                           ve_cases = c(.32, .41, .26, 0),
                           ve_hosps = c(.29, .57, .38, 0),
                           ve_deaths = c(.52, .75, .64, 0))

###################################################################################################
#Regression Model: Simple Model using only age_group and vaccination status as predictors
calibration_data <- covid_hosp %>% filter(weeks_since_july2022 %in% c(0:26))

#Model Calibration (on 6 months of data)
model1 <- glm(num_hosp ~ age_group + boost_2_vax_status_match, 
              data = calibration_data, "quasipoisson")
summary(model1)

#Model Calibration Check (Predicting on 6 months of data that it was calibrated on)
validation_period <- covid_hosp %>% filter(weeks_since_july2022 %in% c(27:52))
pred <- predict(model1, newdata = validation_period, se.fit = T)
validation_period$prediction <- exp(pred$fit)

sum(validation_period$prediction)

###################################################################################################
#Calculating Averted Outcomes

#Creating predictions df
predictions_df <- merge(validation_period %>% group_by(age_group, boost_2_vax_status_match) %>% 
                          summarise(num_hosp = sum(num_hosp), predicted_hosp = sum(prediction)),
                        ve_estimates, by.x = "boost_2_vax_status_match", by.y = "baseline_vacc", all.x = TRUE)

##########################
#Primary Vaccine Strategies Analysis

strat1 <- predictions_df %>% mutate(averted_hosp = ifelse((boost_2_vax_status_match %in% c("Unvaccinated", "Primary Series", "Boosted (1 dose)", "Boosted (2 doses)")),
                                                          predicted_hosp * ve_hosps  , 0) )

strat1$strategy <- rep("strat1", nrow(strat1))
strat2 <- predictions_df %>% mutate(averted_hosp = ifelse((boost_2_vax_status_match %in% c("Primary Series", "Boosted (1 dose)", "Boosted (2 doses)")),
                                                          predicted_hosp * ve_hosps , 0))

strat2$strategy <- rep("strat2", nrow(strat2))
strat3 <- predictions_df %>% mutate(averted_hosp = ifelse((boost_2_vax_status_match %in% c("Unvaccinated")),
                                                          predicted_hosp * ve_hosps ,0))
strat3$strategy <- rep("strat3", nrow(strat3))
strat4 <- predictions_df %>% mutate(averted_hosp = ifelse((boost_2_vax_status_match %in% c("Primary Series")),
                                                          predicted_hosp * ve_hosps , 0))
strat4$strategy <- rep("strat4", nrow(strat3))
strat5 <- predictions_df %>% mutate(averted_hosp = ifelse(((age_group %in% c( "75-84 years", "85+ years")) & (boost_2_vax_status_match != "Unvaccinated")),
                                                          predicted_hosp * ve_hosps ,0))
strat5$strategy <- rep("strat5", nrow(strat3))
strat6 <- predictions_df %>% mutate(averted_hosp = ifelse(((age_group %in% c( "65-74 years", "75-84 years", "85+ years")) & (boost_2_vax_status_match != "Unvaccinated")),
                                                          predicted_hosp * ve_hosps ,0))
strat6$strategy <- rep("strat6", nrow(strat1))
strat7 <- predictions_df %>% mutate(averted_hosp = ifelse(((age_group %in% c( "50-64 years", "65-74 years", "75-84 years", "85+ years")) & (boost_2_vax_status_match != "Unvaccinated")),
                                                          predicted_hosp * ve_hosps, 0))
strat7$strategy <- rep("strat7", nrow(strat1))

combined_strat <- bind_rows(strat1, strat2, strat3, strat4, strat5, strat6, strat7)
##########################
#Additional Vaccine Strategies (Stratified)
strat8 <- predictions_df %>% mutate(averted_hosp = ifelse(((age_group %in% c( "50-64 years", "65-74 years", "75-84 years", "85+ years")) & (boost_2_vax_status_match == "Unvaccinated")),
                                                          predicted_death * ve_deaths, 0))
strat8$strategy <- rep("strat8", nrow(strat8))

strat11 <- predictions_df %>% mutate(averted_death = ifelse(((age_group %in% c( "50-64 years", "65-74 years", "75-84 years", "85+ years")) & (boost_2_vax_status_match == "Primary Series")),
                                                            predicted_death * ve_deaths, 0))
strat11$strategy <- rep("strat11", nrow(strat11))

strat14 <- predictions_df %>% mutate(averted_death = ifelse(((age_group %in% c( "50-64 years", "65-74 years", "75-84 years", "85+ years")) & (boost_2_vax_status_match %in% c("Boosted (1 dose)", "Boosted (2 doses)"))),
                                                            predicted_death * ve_deaths, 0))
strat14$strategy <- rep("strat14", nrow(strat14))


strat9 <- predictions_df %>% mutate(averted_death = ifelse(((age_group %in% c("65-74 years", "75-84 years", "85+ years")) & (boost_2_vax_status_match == "Unvaccinated")),
                                                           predicted_death * ve_deaths, 0))
strat9$strategy <- rep("strat9", nrow(strat9))

strat12 <- predictions_df %>% mutate(averted_death = ifelse(((age_group %in% c( "65-74 years", "75-84 years", "85+ years")) & (boost_2_vax_status_match == "Primary Series")),
                                                            predicted_death * ve_deaths, 0))
strat12$strategy <- rep("strat12", nrow(strat12))

strat15 <- predictions_df %>% mutate(averted_death = ifelse(((age_group %in% c( "65-74 years", "75-84 years", "85+ years")) & (boost_2_vax_status_match %in% c("Boosted (1 dose)", "Boosted (2 doses)"))),
                                                            predicted_death * ve_deaths, 0))
strat15$strategy <- rep("strat15", nrow(strat15))


strat10 <- predictions_df %>% mutate(averted_death = ifelse(((age_group %in% c( "75-84 years", "85+ years")) & (boost_2_vax_status_match == "Unvaccinated")),
                                                            predicted_death * ve_deaths, 0))
strat10$strategy <- rep("strat10", nrow(strat10))

strat13 <- predictions_df %>% mutate(averted_death = ifelse(((age_group %in% c( "75-84 years", "85+ years")) & (boost_2_vax_status_match == "Primary Series")),
                                                            predicted_death * ve_deaths, 0))
strat13$strategy <- rep("strat13", nrow(strat13))

strat16 <- predictions_df %>% mutate(averted_death = ifelse(((age_group %in% c( "75-84 years", "85+ years")) & (boost_2_vax_status_match %in% c("Boosted (1 dose)", "Boosted (2 doses)"))),
                                                            predicted_death * ve_deaths, 0))
strat16$strategy <- rep("strat16", nrow(strat16))


#combined_strat <- bind_rows(strat8, strat9, strat10, strat11, strat12, strat13, strat14, strat15, strat16)

##########################
combined_strat_outcomes <- combined_strat %>% group_by(strategy) %>% summarise(averted_outcomes = ceiling(sum(averted_hosp))) 

averted_outcomes <- melt(combined_strat_outcomes %>% select(strategy, averted_outcomes)) %>% rename(averted_outcomes = value)

merged_data <- merge(averted_outcomes, bivalent_booster_data %>% filter(variable == "bivalent_doses"), by = c("strategy"), all.X = TRUE) %>%
  mutate(predicted_outcomes_baseline = 5650,
         predicted_outcomes_with_vaccine = predicted_outcomes_baseline - averted_outcomes,
         averted_outcome_by_dose = round(averted_outcomes/total_doses * 100000, 2),
         ARR = averted_outcomes/total_doses,
         NNT = ceiling(1/ARR),
         perc_averted = round(averted_outcomes/predicted_outcomes_baseline * 100, 1),
         strategy = factor(strategy, levels = c( "strat1", "strat2", "strat3", "strat4", "strat5", "strat6", "strat7")))

"strat8", "strat9", "strat10", "strat11", "strat12", "strat13", "strat14", "strat15", "strat16")))






###################################################################################################
#Waning VE
predictions_df <- calibration_period %>% group_by(weeks_since_july2022, age_group, boost_2_vax_status_match) %>% summarise(num_cases = sum(num_cases),
                                                                                                                           predicted_cases = sum(prediction)) %>% 
  mutate(ve_cases = case_when((weeks_since_july2022 %in% c(0:4) & boost_2_vax_status_match == 'Unvaccinated') ~ ve$unvax_0_to_4,
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

###################################################################################################
#Prospective Vaccination
predictions_df <- merge(calibration_period %>% group_by(weeks_since_july2022, age_group, boost_2_vax_status_match) %>% summarise(num_cases = sum(num_cases),
                                                                                                                                 predicted_cases = sum(prediction)) %>% 
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
                                                           predicted_cases * ve_cases* vax_coverage, 0))

strat2$strategy <- rep("strat2", nrow(strat2))
strat3 <- predictions_df %>% mutate(averted_cases = ifelse((boost_2_vax_status_match %in% c("Unvaccinated")),
                                                           predicted_cases * ve_cases* vax_coverage,0))
strat3$strategy <- rep("strat3", nrow(strat3))
strat4 <- predictions_df %>% mutate(averted_cases = ifelse((boost_2_vax_status_match %in% c("Primary Series")),
                                                           predicted_cases * ve_cases* vax_coverage, 0))
strat4$strategy <- rep("strat4", nrow(strat3))
strat5 <- predictions_df %>% mutate(averted_cases = ifelse(((age_group %in% c( "75-84 years", "85+ years")) & (boost_2_vax_status_match != "Unvaccinated")),
                                                           predicted_cases * ve_cases* vax_coverage,0))
strat5$strategy <- rep("strat5", nrow(strat3))
strat6 <- predictions_df %>% mutate(averted_cases = ifelse(((age_group %in% c( "65-74 years", "75-84 years", "85+ years")) & (boost_2_vax_status_match != "Unvaccinated")),
                                                           predicted_cases * ve_cases* vax_coverage,0))
strat6$strategy <- rep("strat6", nrow(strat1))
strat7 <- predictions_df %>% mutate(averted_cases = ifelse(((age_group %in% c( "50-64 years", "65-74 years", "75-84 years", "85+ years")) & (boost_2_vax_status_match != "Unvaccinated")),
                                                           predicted_cases * ve_cases* vax_coverage, 0))
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
         NNT = ceiling(1/ARR),
         perc_averted = round(averted_outcomes/predicted_outcomes_baseline * 100, 1),
         strategy = factor(strategy, levels = c("strat1", "strat2", "strat3", "strat4", "strat5", "strat6", "strat7"))) %>%
  arrange(strategy)



