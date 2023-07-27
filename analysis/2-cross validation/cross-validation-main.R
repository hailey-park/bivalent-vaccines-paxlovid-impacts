###################################################################################################
#Title: Cross Validation
#Author: Hailey Park
#Date: January 30, 2023
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
adjustment_factor_cases <- read.csv("outcome data clean/adjustment_factor_matrix_cases_crossval.csv")[,-1]
adjustment_factor_hosp <- read.csv("outcome data clean/adjustment_factor_matrix_hosp_crossval.csv")[,-1]
adjustment_factor_death <- read.csv("outcome data clean/adjustment_factor_matrix_death_crossval.csv")[,-1]
covid_case_data_crossval <- read.csv("outcome data clean/covid_data_clean_cases_crossval.csv")[,-1]
covid_hosp_data_crossval <- read.csv("outcome data clean/covid_data_clean_hosp_crossval.csv")[,-1]
covid_death_data_crossval <- read.csv("outcome data clean/covid_data_clean_death_crossval.csv")[,-1]


#Make dataframe for month to week conversion (later used for merging adj factor df)
month_weeks <- covid_death_data_crossval %>% 
  mutate(weeks_since_april2022 = as.integer(floor(difftime(strptime(dtepisode, format="%Y-%m-%d"), 
                                                          strptime("2022-04-23", format="%Y-%m-%d"), units = 'weeks'))),
         month = as.factor(floor_date(as.Date(dtepisode), "month"))) %>%
  group_by(weeks_since_april2022, month) %>% summarise(month = names(table(month))[which.max(table(month))]) %>% group_by(weeks_since_april2022) %>% summarise(month = names(table(month))[which.max(table(month))])

#Make dataframes for each outcome (cases, hosp, deaths)
covid_cases <- merge(merge(covid_case_data_crossval %>% 
                             mutate(weeks_since_april2022 = as.integer(floor(difftime(strptime(dtepisode, format="%Y-%m-%d"), 
                                                                                     strptime("2022-04-23", format="%Y-%m-%d"), units = 'weeks')))) %>%
                             group_by(weeks_since_april2022, age_group, boost_2_vax_status_match, .drop = FALSE) %>% summarise(num_cases = n()),
                           month_weeks, by = "weeks_since_april2022", all.x = TRUE),
                     adjustment_factor_cases, by = c("month", "age_group", "boost_2_vax_status_match"), all.x = TRUE) %>%
  mutate(num_cases_adj = ceiling(num_cases * adj_factor)) %>% select(weeks_since_april2022, boost_2_vax_status_match, age_group, num_cases, adj_factor, num_cases_adj)



covid_hosp <- merge(merge(covid_hosp_data_crossval %>% 
                            mutate(weeks_since_april2022 = as.integer(floor(difftime(strptime(dtepisode, format="%Y-%m-%d"), 
                                                                                    strptime("2022-04-23", format="%Y-%m-%d"), units = 'weeks')))) %>%
                            group_by(weeks_since_april2022, age_group, boost_2_vax_status_match, .drop = FALSE) %>% summarise(num_hosp = n()),
                          month_weeks, by = "weeks_since_april2022", all.x = TRUE),
                    adjustment_factor_hosp, by = c("month", "age_group", "boost_2_vax_status_match"), all.x = TRUE) %>%
  mutate(num_hosp_adj = ceiling(num_hosp * adj_factor * (1/0.75))) %>% select(weeks_since_april2022, boost_2_vax_status_match, age_group, num_hosp, adj_factor, num_hosp_adj)


covid_death <- merge(merge(covid_death_data_crossval %>% 
                             mutate(weeks_since_april2022 = as.integer(floor(difftime(strptime(dtepisode, format="%Y-%m-%d"), 
                                                                                     strptime("2022-04-23", format="%Y-%m-%d"), units = 'weeks')))) %>%
                             group_by(weeks_since_april2022, age_group, boost_2_vax_status_match, .drop = FALSE) %>% summarise(num_death = n()),
                           month_weeks, by = "weeks_since_april2022", all.x = TRUE),
                     adjustment_factor_death, by = c("month", "age_group", "boost_2_vax_status_match"), all.x = TRUE) %>%
  mutate(num_death_adj = ceiling(num_death * adj_factor)) %>% select(weeks_since_april2022, boost_2_vax_status_match, age_group, num_death, adj_factor, num_death_adj)

###################################################################################################
#Regression Model: Simple Model using only age_group and vaccination status as predictors
#NOTE: We are running the regression model on each outcome (cases, hospitalizations, deaths) separately, so make sure to swap out dependent variable ('num_cases_adj', 'num_hosp_adj', 'num_death_adj'), 
#      and outcome dataframes ('covid_case', 'covid_hosp', 'covid_death').

calibration_data <- covid_cases %>% filter(weeks_since_april2022 %in% c(0:26))

#Model Calibration (on 6 months of data)
model1 <- glm(num_cases_adj ~ age_group + boost_2_vax_status_match, 
              data = calibration_data, "quasipoisson")
summary(model1)

#Model Calibration Check (Predicting on 6 months of data that it was calibrated on)
validation_period <- covid_cases %>% filter(weeks_since_april2022 %in% c(27:39))
pred <- predict(model1, newdata = validation_period, se.fit = T)
validation_period$prediction <- exp(pred$fit)

sum(validation_period$prediction)
sum((covid_cases%>% filter(weeks_since_april2022 %in% c(27:39)))$num_cases_adj, na.rm = TRUE)
###################################################################################################
#Calculating Predicted Outcomes by Strategy

#Creating predictions df
predictions_df <- validation_period %>% group_by(age_group, boost_2_vax_status_match) %>% 
                                summarise(num_cases = sum(num_cases_adj), predicted_cases = sum(prediction))

#Primary Vaccine Strategies Analysis

strat1 <- predictions_df %>% mutate(predicted_cases = ifelse((boost_2_vax_status_match %in% c("Unvaccinated", "Primary Series", "Boosted (1 dose)", "Boosted (2 doses)")),
                                                           predicted_cases, 0),
                                    num_cases = ifelse((boost_2_vax_status_match %in% c("Unvaccinated", "Primary Series", "Boosted (1 dose)", "Boosted (2 doses)")),
                                                       num_cases, 0))

strat1$strategy <- rep("strat1", nrow(strat1))
strat2 <- predictions_df %>% mutate(predicted_cases = ifelse((boost_2_vax_status_match %in% c("Primary Series", "Boosted (1 dose)", "Boosted (2 doses)")),
                                                           predicted_cases, 0),
                                    num_cases = ifelse((boost_2_vax_status_match %in% c("Primary Series", "Boosted (1 dose)", "Boosted (2 doses)")),
                                                       num_cases, 0))

strat2$strategy <- rep("strat2", nrow(strat2))
strat3 <- predictions_df %>% mutate(predicted_cases = ifelse((boost_2_vax_status_match %in% c("Unvaccinated")),
                                                           predicted_cases,0),
                                    num_cases = ifelse((boost_2_vax_status_match %in% c("Unvaccinated")),
                                                       num_cases,0))

strat3$strategy <- rep("strat3", nrow(strat3))
strat4 <- predictions_df %>% mutate(predicted_cases = ifelse((boost_2_vax_status_match %in% c("Primary Series")),
                                                           predicted_cases, 0),
                                    num_cases = ifelse((boost_2_vax_status_match %in% c("Primary Series")),
                                                                    num_cases, 0))

strat4$strategy <- rep("strat4", nrow(strat3))
strat5 <- predictions_df %>% mutate(predicted_cases = ifelse(((age_group %in% c( "75-84 years", "85+ years")) & (boost_2_vax_status_match != "Unvaccinated")),
                                                           predicted_cases,0),
                                    num_cases =  ifelse(((age_group %in% c( "75-84 years", "85+ years")) & (boost_2_vax_status_match != "Unvaccinated")),
                                                                     num_cases,0)) 
strat5$strategy <- rep("strat5", nrow(strat3))
strat6 <- predictions_df %>% mutate(predicted_cases = ifelse(((age_group %in% c( "65-74 years", "75-84 years", "85+ years")) & (boost_2_vax_status_match != "Unvaccinated")),
                                                           predicted_cases,0),
                                    num_cases = ifelse(((age_group %in% c( "65-74 years", "75-84 years", "85+ years")) & (boost_2_vax_status_match != "Unvaccinated")),
                                                                   num_cases,0))
strat6$strategy <- rep("strat6", nrow(strat1))
strat7 <- predictions_df %>% mutate(predicted_cases = ifelse(((age_group %in% c( "50-64 years", "65-74 years", "75-84 years", "85+ years")) & (boost_2_vax_status_match != "Unvaccinated")),
                                                           predicted_cases, 0),
                                    num_cases = ifelse(((age_group %in% c( "50-64 years", "65-74 years", "75-84 years", "85+ years")) & (boost_2_vax_status_match != "Unvaccinated")),
                                                                    num_cases, 0))
strat7$strategy <- rep("strat7", nrow(strat1))

combined_strat <- bind_rows(strat1, strat2, strat3, strat4, strat5, strat6, strat7)


combined_strat_outcomes <- combined_strat %>% group_by(strategy) %>% summarise(predicted_cases = round(sum(predicted_cases)),
                                                                               num_cases = round(sum(num_cases))
) 

