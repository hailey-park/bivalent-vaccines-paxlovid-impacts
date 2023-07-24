###################################################################################################
#Title: Bivalent Booster Regression Analysis
#Author: Hailey Park
#Date: July 3, 2023
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
covid_case_data <- read.csv("outcome data clean/covid_data_clean_cases_final.csv")[,-1]
covid_hosp_data <- read.csv("outcome data clean/covid_data_clean_hosp_final.csv")[,-1]
covid_death_data <- read.csv("outcome data clean/covid_data_clean_death_final.csv")[,-1]
bivalent_booster_data <- read.csv("outcome data clean/bivalent_booster_doses_clean_70maxuptake_v2.csv")[,-1]
adjustment_factor_cases <- read.csv("outcome data clean/adjustment_factor_matrix_cases-sensitivity-oldest age bivalent coverage.csv")[,-1]
adjustment_factor_hosp <- read.csv("outcome data clean/adjustment_factor_matrix_hosp-sensitivity-oldest age bivalent coverage.csv")[,-1]
adjustment_factor_death <- read.csv("outcome data clean/adjustment_factor_matrix_death-sensitivity-oldest age bivalent coverage.csv")[,-1]
bivalent_coverage <- read.csv("outcome data clean/bivalent_coverage_age_vax-sensitivity-oldest age bivalent coverage.csv")[,-1]

  
#Make dataframe for month to week conversion (later used for merging adj factor df)
month_weeks <- covid_death_data %>% 
  mutate(weeks_since_july2022 = as.integer(floor(difftime(strptime(dtepisode, format="%Y-%m-%d"),strptime("2022-07-23", format="%Y-%m-%d"), units = 'weeks'))),
         month = as.factor(floor_date(as.Date(dtepisode), "month"))) %>%
  group_by(weeks_since_july2022, month) %>% summarise(month = names(table(month))[which.max(table(month))]) %>% group_by(weeks_since_july2022) %>% summarise(month = names(table(month))[which.max(table(month))])

#Make dataframes for each outcome (cases, hosp, deaths)
covid_cases <- merge(merge(covid_case_data %>% mutate(weeks_since_july2022 = as.integer(floor(difftime(strptime(dtepisode, format="%Y-%m-%d"), 
                                                                                                       strptime("2022-07-23", format="%Y-%m-%d"), units = 'weeks')))) %>%
                             group_by(weeks_since_july2022, age_group, boost_2_vax_status_match, .drop = FALSE) %>% summarise(num_cases = n()),
                           month_weeks, by = "weeks_since_july2022", all.x = TRUE),
                     adjustment_factor_cases, by = c("month", "age_group", "boost_2_vax_status_match"), all.x = TRUE) %>%
  mutate(num_cases_adj = ceiling(num_cases * adj_factor)) %>% select(weeks_since_july2022, boost_2_vax_status_match, age_group, num_cases, adj_factor, num_cases_adj)


covid_hosp <- merge(merge(covid_hosp_data %>% mutate(weeks_since_july2022 = as.integer(floor(difftime(strptime(dtepisode, format="%Y-%m-%d"), 
                                                                                                      strptime("2022-07-23", format="%Y-%m-%d"), units = 'weeks')))) %>%
                            group_by(weeks_since_july2022, age_group, boost_2_vax_status_match, .drop = FALSE) %>% summarise(num_hosp = n()),
                          month_weeks, by = "weeks_since_july2022", all.x = TRUE),
                    adjustment_factor_hosp, by = c("month", "age_group", "boost_2_vax_status_match"), all.x = TRUE) %>%
  mutate(num_hosp_adj = ceiling(num_hosp * adj_factor * (1/0.75))) %>% select(weeks_since_july2022, boost_2_vax_status_match, age_group, num_hosp, adj_factor, num_hosp_adj)

  
covid_death <- merge(merge(covid_death_data %>% mutate(weeks_since_july2022 = as.integer(floor(difftime(strptime(dtepisode, format="%Y-%m-%d"), strptime("2022-07-23", format="%Y-%m-%d"), units = 'weeks')))) %>%
                             group_by(weeks_since_july2022, age_group, boost_2_vax_status_match, .drop = FALSE) %>% summarise(num_death = n()),
                           month_weeks, by = "weeks_since_july2022", all.x = TRUE),
                     adjustment_factor_death, by = c("month", "age_group", "boost_2_vax_status_match"), all.x = TRUE) %>%
  mutate(num_death_adj = ceiling(num_death * adj_factor)) %>% select(weeks_since_july2022, boost_2_vax_status_match, age_group, num_death, adj_factor, num_death_adj)

  
#Create dataframes for VE estimates for averting outcomes
ve_estimates <- data.frame(baseline_vacc = c("Unvaccinated", "Primary Series", "Boosted (1 dose)","Boosted (2 doses)"),
                           #ve_cases = c(.32, .41, .26, .40),
                           ve_hosps = c(.29, .57, .38, .56),
                           ve_deaths = c(.57, .75, .62, .63))

###################################################################################################
#Regression Model: Simple Model using only age_group and vaccination status as predictors
calibration_data <- covid_cases %>% filter(weeks_since_july2022 %in% c(0:26))

#Model Calibration (on 6 months of data)
model1 <- glm(num_cases_adj ~ age_group + boost_2_vax_status_match,
              data = calibration_data, "quasipoisson")
summary(model1)


#Model Calibration Check (Predicting on 6 months of data that it was calibrated on)
validation_period <- covid_cases %>% filter(weeks_since_july2022 %in% c(0:26))
pred <- predict(model1, newdata = validation_period, se.fit = T)
validation_period$prediction <- exp(pred$fit)

sum(validation_period$prediction)
sum(covid_cases$num_cases_adj)

dat1 <- add_pi(validation_period, model1, names = c("lpb", "upb"), alpha = 0.1, nsims = 10000)
dat2 <- add_ci(validation_period, model1, names = c("lpb", "upb"))
                
                  
sum(dat2$pred)
sum(dat2$lpb)
sum(dat2$upb)

#For hospitalizations/deaths
predictions_df <- merge(merge(validation_period %>% group_by(age_group, boost_2_vax_status_match) %>% 
                                summarise(num_hosp = sum(num_hosp_adj), predicted_hosp = sum(prediction)),
                              ve_estimates, by.x = "boost_2_vax_status_match", by.y = "baseline_vacc", all.x = TRUE),
                        bivalent_coverage, by = c("boost_2_vax_status_match", "age_group"), all.x = TRUE) %>%
  select(-num_hosp) %>% mutate(predicted_hosp = round(predicted_hosp))

#For cases (add waning VE estimates)
predictions_df <- merge(validation_period %>% group_by(weeks_since_july2022, age_group, boost_2_vax_status_match) %>% 
                          summarise(num_cases = sum(num_cases_adj),
                                    predicted_cases = sum(prediction)) %>% 
                          mutate(ve_cases = case_when((weeks_since_july2022 %in% c(0:4) & boost_2_vax_status_match == 'Unvaccinated') ~ 0.45,
                                                      (weeks_since_july2022 %in% c(5:26) & boost_2_vax_status_match == 'Unvaccinated') ~ 0.32,
                                                      (weeks_since_july2022 %in% c(0:4) & boost_2_vax_status_match == 'Primary Series') ~ 0.56,
                                                      (weeks_since_july2022 %in% c(5:9) & boost_2_vax_status_match == 'Primary Series') ~ 0.43,
                                                      (weeks_since_july2022 %in% c(10:13) & boost_2_vax_status_match == 'Primary Series') ~ 0.26,
                                                      (weeks_since_july2022 %in% c(14:17) & boost_2_vax_status_match == 'Primary Series') ~ 0.20,
                                                      (weeks_since_july2022 %in% c(18:26) & boost_2_vax_status_match == 'Primary Series') ~ 0.19,
                                                      (weeks_since_july2022 %in% c(0:5) & boost_2_vax_status_match == 'Boosted (1 dose)') ~ 0.39,
                                                      (weeks_since_july2022 %in% c(6:15) & boost_2_vax_status_match == 'Boosted (1 dose)') ~ 0.32,
                                                      (weeks_since_july2022 %in% c(16:26) & boost_2_vax_status_match == 'Boosted (1 dose)') ~ 0,
                                                      (weeks_since_july2022 %in% c(0:5) & boost_2_vax_status_match == 'Boosted (2 doses)') ~ 0.39,
                                                      (weeks_since_july2022 %in% c(6:15) & boost_2_vax_status_match == 'Boosted (2 doses)') ~ 0.32,
                                                      (weeks_since_july2022 %in% c(16:26) & boost_2_vax_status_match == 'Boosted (2 doses)') ~ 0,
                                                      TRUE ~ 0)),
                        bivalent_coverage, by = c("boost_2_vax_status_match", "age_group"), all.x = TRUE) %>%
  mutate(predicted_cases = round(predicted_cases))

##########################
#Primary Vaccine Strategies Analysis
                                  
strat1 <- predictions_df %>% mutate(averted_diff = ifelse((boost_2_vax_status_match %in% c("Unvaccinated", "Primary Series", "Boosted (1 dose)", "Boosted (2 doses)")),
                                                          predicted_cases * ve_cases * (0.7 - bivalent_coverage), 0),
                                    total_predicted = predicted_cases * (1 - (bivalent_coverage * ve_cases)))

strat1$strategy <- rep("strat1", nrow(strat1))

strat2 <- predictions_df %>% mutate(averted_diff = ifelse((boost_2_vax_status_match %in% c("Primary Series", "Boosted (1 dose)", "Boosted (2 doses)")),
                                                          predicted_cases * ve_cases * (0.7 - bivalent_coverage), 0),
                                    total_predicted = predicted_cases * (1 - (bivalent_coverage * ve_cases)))
strat2$strategy <- rep("strat2", nrow(strat2))

strat3 <- predictions_df %>% mutate(averted_diff = ifelse((boost_2_vax_status_match %in% c("Unvaccinated")),
                                                          predicted_cases * ve_cases * (0.7 - bivalent_coverage), 0),
                                    total_predicted = predicted_cases * (1 - (bivalent_coverage * ve_cases)))
strat3$strategy <- rep("strat3", nrow(strat3))

strat4 <- predictions_df %>% mutate(averted_diff = ifelse((boost_2_vax_status_match %in% c("Primary Series")),
                                                          predicted_cases * ve_cases * (0.7 - bivalent_coverage), 0),
                                    total_predicted = predicted_cases * (1 - (bivalent_coverage * ve_cases)))
strat4$strategy <- rep("strat4", nrow(strat3))

strat5 <- predictions_df %>% mutate(averted_diff = ifelse(((age_group %in% c( "75-84 years", "85+ years")) & (boost_2_vax_status_match != "Unvaccinated")),
                                                          predicted_cases * ve_cases * (0.7 - bivalent_coverage), 0),
                                    total_predicted = predicted_cases * (1 - (bivalent_coverage * ve_cases))) 
strat5$strategy <- rep("strat5", nrow(strat3))

strat6 <- predictions_df %>% mutate(averted_diff = ifelse(((age_group %in% c( "65-74 years", "75-84 years", "85+ years")) & (boost_2_vax_status_match != "Unvaccinated")),
                                                          predicted_cases * ve_cases * (0.7 - bivalent_coverage), 0),
                                    total_predicted = predicted_cases * (1 - (bivalent_coverage * ve_cases)))
strat6$strategy <- rep("strat6", nrow(strat1))

strat7 <- predictions_df %>% mutate(averted_diff = ifelse(((age_group %in% c("50-64 years", "65-74 years", "75-84 years", "85+ years")) & (boost_2_vax_status_match != "Unvaccinated")),
                                                          predicted_cases * ve_cases * (0.7 - bivalent_coverage), 0),
                                    total_predicted = predicted_cases * (1 - (bivalent_coverage * ve_cases)))
strat7$strategy <- rep("strat7", nrow(strat1))

                                
##########################

combined_strat <- bind_rows(strat1, strat2, strat3, strat4, strat5, strat6, strat7)


combined_strat_outcomes <- combined_strat %>% group_by(strategy) %>% summarise(averted_outcomes_diff = round(sum(averted_diff)),
                                                                               total_predicted_outcomes = round(sum(total_predicted))) 


merged_data <- merge(combined_strat_outcomes, bivalent_booster_data, by = c("strategy"), all.X = TRUE) %>%
  mutate(ARR = averted_outcomes_diff/(total_doses),
         NNT = ceiling(1/ARR),
         perc_averted = round(averted_outcomes_diff/total_predicted_outcomes * 100, 1),
         strategy = factor(strategy, levels = c("strat1", "strat2", "strat3", "strat4", "strat5", "strat6", "strat7")))%>%
  select(-c(variable))



