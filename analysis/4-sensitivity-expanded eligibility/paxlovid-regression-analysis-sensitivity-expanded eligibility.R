###################################################################################################
#Title: Nirmatrelvir-Ritonavir Regression Analysis
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
covid_case_data <- read.csv("outcome data clean/covid_data_clean_cases_final.csv")[,-1]
covid_hosp_data <- read.csv("outcome data clean/covid_data_clean_hosp_final.csv")[,-1]
covid_death_data <- read.csv("outcome data clean/covid_data_clean_death_final.csv")[,-1]
adjustment_factor_cases <- read.csv("outcome data clean/adjustment_factor_matrix_cases.csv")[,-1]
adjustment_factor_hosp <- read.csv("outcome data clean/adjustment_factor_matrix_hosp.csv")[,-1]
adjustment_factor_death <- read.csv("outcome data clean/adjustment_factor_matrix_death.csv")[,-1]

#Make dataframe for month to week conversion (later used for merging adj factor df)
month_weeks <- covid_death_data %>% 
  mutate(weeks_since_july2022 = as.integer(floor(difftime(strptime(dtepisode, format="%Y-%m-%d"), 
                                                          strptime("2022-07-23", format="%Y-%m-%d"), units = 'weeks'))),
         month = as.factor(floor_date(as.Date(dtepisode), "month"))) %>%
  group_by(weeks_since_july2022, month) %>% summarise(month = names(table(month))[which.max(table(month))]) %>% group_by(weeks_since_july2022) %>% summarise(month = names(table(month))[which.max(table(month))])

#Make dataframes for each outcome (cases, hosp, deaths)
covid_cases <- merge(merge(covid_case_data %>% 
                             mutate(weeks_since_july2022 = as.integer(floor(difftime(strptime(dtepisode, format="%Y-%m-%d"), 
                                                                                     strptime("2022-07-23", format="%Y-%m-%d"), units = 'weeks')))) %>%
                             group_by(weeks_since_july2022, age_group, boost_2_vax_status_match, .drop = FALSE) %>% summarise(num_cases = n()),
                           month_weeks, by = "weeks_since_july2022", all.x = TRUE),
                     adjustment_factor_cases, by = c("month", "age_group", "boost_2_vax_status_match"), all.x = TRUE) %>%
  mutate(num_cases_adj = ceiling(num_cases * adj_factor)) %>% select(weeks_since_july2022, boost_2_vax_status_match, age_group, num_cases, adj_factor, num_cases_adj)



covid_hosp <- merge(merge(covid_hosp_data %>% 
                            mutate(weeks_since_july2022 = as.integer(floor(difftime(strptime(dtepisode, format="%Y-%m-%d"), 
                                                                                    strptime("2022-07-23", format="%Y-%m-%d"), units = 'weeks')))) %>%
                            group_by(weeks_since_july2022, age_group, boost_2_vax_status_match, .drop = FALSE) %>% summarise(num_hosp = n()),
                          month_weeks, by = "weeks_since_july2022", all.x = TRUE),
                    adjustment_factor_hosp, by = c("month", "age_group", "boost_2_vax_status_match"), all.x = TRUE) %>%
  mutate(num_hosp_adj = ceiling(num_hosp * adj_factor * (1/0.75))) %>% select(weeks_since_july2022, boost_2_vax_status_match, age_group, num_hosp, adj_factor, num_hosp_adj)


covid_death <- merge(merge(covid_death_data %>% 
                             mutate(weeks_since_july2022 = as.integer(floor(difftime(strptime(dtepisode, format="%Y-%m-%d"), 
                                                                                     strptime("2022-07-23", format="%Y-%m-%d"), units = 'weeks')))) %>%
                             group_by(weeks_since_july2022, age_group, boost_2_vax_status_match, .drop = FALSE) %>% summarise(num_death = n()),
                           month_weeks, by = "weeks_since_july2022", all.x = TRUE),
                     adjustment_factor_death, by = c("month", "age_group", "boost_2_vax_status_match"), all.x = TRUE) %>%
  mutate(num_death_adj = ceiling(num_death * adj_factor)) %>% select(weeks_since_july2022, boost_2_vax_status_match, age_group, num_death, adj_factor, num_death_adj)

#Merge the case outcome df with the severe outcome dfs
pax_death <- merge(covid_cases, covid_death, by = c("weeks_since_july2022", "age_group", "boost_2_vax_status_match"),
                   all.x = TRUE) %>% replace(is.na(.), 0)

pax_hosp <- merge(covid_cases, covid_hosp, by = c("weeks_since_july2022", "age_group", "boost_2_vax_status_match"), 
                  all.x = TRUE)


pe_estimates <- data.frame(baseline_vacc = c("Unvaccinated", "Vaccinated"),
                               pe_hosps = c(.89, .40),
                               pe_deaths = c(.89, .71))



inspection <- covid_cases %>% group_by(age_group, boost_2_vax_status_match) %>% summarise(total_death = sum(num_cases_adj))

###################################################################################################
#Regression Model: Simple Model using only age_group and vaccination status as predictors
#NOTE: We are running the regression model on each outcome (cases, hospitalizations, deaths) separately, so make sure to swap out dependent variable ('num_cases_adj', 'num_hosp_adj', 'num_death_adj'), 
#      the vaccine effectiveness names ('pe_hosps', 'pe_deaths') and outcome dataframes ('pax_hosp', 'pax_death').


#Model Calibration (on 6 months of data)
model1 <- glm(num_cases_adj ~ age_group + boost_2_vax_status_match, 
              data = pax_death %>% filter(weeks_since_july2022 %in% c(0:26)), "quasipoisson")
summary(model1)

#Model Calibration Check (Predicting on 6 months of data that it was calibrated on)
calibration_period <- pax_death %>% filter(weeks_since_july2022 %in% c(0:26))
pred <- predict(model1, newdata = calibration_period, se.fit = T)
calibration_period$prediction_cases <- exp(pred$fit)

sum(calibration_period$num_death_adj)
sum(calibration_period$prediction_death)
###################################################################################################
#Calculating Averted Outcomes

#Creating predictions df
predictions_pax_df <- merge(calibration_period %>% mutate(vacc_status = ifelse(boost_2_vax_status_match %in% c("Primary Series", "Boosted (1 dose)", "Boosted (2 doses)"),
                                                                               "Vaccinated",
                                                                               "Unvaccinated")) %>% group_by(age_group, boost_2_vax_status_match, vacc_status) %>% 
                              summarise(num_cases = sum(num_cases_adj), 
                                        predicted_cases = sum(prediction_cases), 
                                        num_death = sum(num_death_adj),
                                        predicted_death = sum(prediction_death)),
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


#Strategies with Expanded Eligibility (Strategies 1-2)

#Strategy 1: Giving Paxlovid to those 18+ years (everyone)
  #Assumption: setting maximum paxlovid coverage to be 70% for everyone
strat1 <- predictions_pax_df %>% mutate(averted_death_diff = case_when((age_group %in% c("18-49 years", "50-64 years", "65-74 years", "75-84 years", "85+ years"))  ~ predicted_death * pe_deaths * (0.7 - baseline_pax_uptake),
                                                                  TRUE ~ 0),
                                        target_group_cases_diff = case_when((age_group %in% c("18-49 years", "50-64 years", "65-74 years", "75-84 years", "85+ years")) ~ predicted_cases * (0.7 - baseline_pax_uptake),
                                                                       TRUE ~ 0),
                                        total_predicted = predicted_death * (1 - (baseline_pax_uptake * pe_deaths)))
strat1$strategy <- rep("strat1", nrow(strat1))

#Strategy 2: Giving Paxlovid to those 50+ years (everyone)
  #Assumption: setting maximum paxlovid coverage to be 70% for everyone
strat2 <- predictions_pax_df %>% mutate(averted_death_diff = case_when((age_group %in% c("50-64 years", "65-74 years", "75-84 years", "85+ years"))  ~ predicted_death * pe_deaths * (0.7 - baseline_pax_uptake),
                                                                  TRUE ~ 0),
                                        target_group_cases_diff = case_when((age_group %in% c("50-64 years", "65-74 years", "75-84 years", "85+ years")) ~ predicted_cases * (0.7 - baseline_pax_uptake),
                                                                       TRUE ~ 0),
                                        total_predicted = predicted_death * (1 - (baseline_pax_uptake * pe_deaths)))
strat2$strategy <- rep("strat2", nrow(strat2))

combined_strat <- bind_rows(strat1, strat2)


combined_strat_outcomes <- combined_strat %>% group_by(strategy) %>% summarise(averted_outcomes_diff = ceiling(sum(averted_death_diff)),
                                                                               total_death = ceiling(sum(total_predicted)),
                                                                               target_group_cases_diff = ceiling(sum(target_group_cases_diff))) %>%
  mutate(averted_outcomes_cascade = ceiling(averted_outcomes_diff * 0.8 * 0.895),
         treated_cascade = ceiling(target_group_cases_diff * 0.8 * 0.895),
         NNT = ceiling(treated_cascade/averted_outcomes_cascade),
         perc_averted = round(averted_outcomes_cascade/ total_death * 100, 1),
         strategy = factor(strategy, levels = c("strat1", "strat2")))
         

