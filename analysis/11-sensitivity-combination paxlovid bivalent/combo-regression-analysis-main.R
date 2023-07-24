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
bivalent_coverage <- read.csv("outcome data clean/bivalent_coverage_age_vax.csv")[,-1]


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


pax_death <- merge(covid_cases, covid_death, by = c("weeks_since_july2022", "age_group", "boost_2_vax_status_match"),
                   all.x = TRUE) %>% replace(is.na(.), 0)

pax_hosp <- merge(covid_cases, covid_hosp, by = c("weeks_since_july2022", "age_group", "boost_2_vax_status_match"), 
                  all.x = TRUE)


pe_estimates <- data.frame(baseline_vacc = c("Unvaccinated", "Vaccinated"),
                               pe_hosps = c(.89, .40),
                               pe_deaths = c(.89, .71))

ve_estimates <- data.frame(baseline_vacc = c("Unvaccinated", "Primary Series", "Boosted (1 dose)","Boosted (2 doses)"),
                           ve_hosps = c(.29, .57, .38, .56),
                           ve_deaths = c(.57, .75, .62, .63))

###################################################################################################
#Regression Model: Simple Model using only age_group and vaccination status as predictors

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


inspection <- covid_cases %>% group_by(age_group) %>% summarise(total_cases = sum(num_cases_adj))
###################################################################################################
#Calculating Averted Outcomes

#Creating predictions df
predictions_pax_df <- merge(merge(merge(calibration_period %>% mutate(vacc_status = ifelse(boost_2_vax_status_match %in% c("Primary Series", "Boosted (1 dose)", "Boosted (2 doses)"),
                                                                               "Vaccinated",
                                                                               "Unvaccinated")) %>% group_by(age_group, boost_2_vax_status_match, vacc_status) %>% 
                              summarise(num_cases = sum(num_cases_adj), 
                                        predicted_cases = sum(prediction_cases), 
                                        num_death = sum(num_death_adj),
                                        predicted_death = sum(prediction_death)),
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


#######################################
#Primary Strategies

#Strategy 3: Giving Paxlovid to those 18+ years (unvaccinated and comorbidities)
strat3 <- predictions_pax_df %>% mutate(averted_death_diff_vax = case_when((age_group %in% c("18-49 years","50-64 years", "65-74 years", "75-84 years", "85+ years")) ~ predicted_death * ve_deaths * (0.7 - bivalent_coverage),
                                                                          TRUE ~ 0), 
                                        averted_death_diff_pax = case_when((age_group %in% c("18-49 years","50-64 years", "65-74 years", "75-84 years", "85+ years")) ~ (predicted_death - averted_death_diff_vax) * pe_deaths * (max_pax_uptake - baseline_pax_uptake),
                                                                 TRUE ~ 0),
                                        total_predicted = predicted_death * (1 - (baseline_pax_uptake * pe_deaths) - (bivalent_coverage * ve_deaths)))
                                        

strat3$strategy <- rep("strat3", nrow(strat3))

#Strategy 4: Giving Paxlovid to those 50+ years (unvaccinated and comorbidities)
strat4 <- predictions_pax_df %>% mutate(averted_death_diff_vax = case_when((age_group %in% c("50-64 years", "65-74 years", "75-84 years", "85+ years")) ~ predicted_death * ve_deaths * (0.7 - bivalent_coverage),
                                                                          TRUE ~ 0),
                                        averted_death_diff_pax = case_when((age_group %in% c("50-64 years", "65-74 years", "75-84 years", "85+ years")) ~ (predicted_death - averted_death_diff_vax) * pe_deaths * (max_pax_uptake - baseline_pax_uptake),
                                                                           TRUE ~ 0),
                                        total_predicted = predicted_death * (1 - (baseline_pax_uptake * pe_deaths) - (bivalent_coverage * ve_deaths)))  

strat4$strategy <- rep("strat4", nrow(strat4))


#Strategy 5: Giving Paxlovid to the those 65+ years
strat5 <- predictions_pax_df %>% mutate(averted_death_diff_vax = case_when((age_group %in% c("65-74 years", "75-84 years", "85+ years")) ~ predicted_death * ve_deaths * (0.7 - bivalent_coverage),
                                                                          TRUE ~ 0),
                                        averted_death_diff_pax = ifelse((age_group %in% c("65-74 years", "75-84 years", "85+ years")),
                                                                       (predicted_death - averted_death_diff_vax) * pe_deaths * (max_pax_uptake - baseline_pax_uptake), 0),
                                        total_predicted = predicted_death * (1 - (baseline_pax_uptake * pe_deaths) - (bivalent_coverage * ve_deaths)))

strat5$strategy <- rep("strat5", nrow(strat5))

#Strategy 6: Giving Paxlovid to the those 75+ years
strat6 <- predictions_pax_df %>% mutate(averted_death_diff_vax = case_when((age_group %in% c("75-84 years", "85+ years")) ~ predicted_death * ve_deaths * (0.7 - bivalent_coverage),
                                                                          TRUE ~ 0),
                                        averted_death_diff_pax = ifelse((age_group %in% c("75-84 years", "85+ years")),
                                                                       (predicted_death - averted_death_diff_vax) * pe_deaths * (max_pax_uptake - baseline_pax_uptake), 0),
                                        total_predicted = predicted_death * (1 - (baseline_pax_uptake * pe_deaths) - (bivalent_coverage * ve_deaths)))

strat6$strategy <- rep("strat6", nrow(strat6))  

combined_strat <- bind_rows(strat3, strat4, strat5, strat6)


#######################################
combined_strat_outcomes <- combined_strat %>% group_by(strategy) %>% summarise(averted_outcomes_diff_vax = ceiling(sum(averted_death_diff_vax)),
                                                                               averted_outcomes_diff_pax = ceiling(sum(averted_death_diff_pax)),
                                                                               total_death = ceiling(sum(total_predicted))) %>%
  mutate(averted_outcomes_cascade_pax = ceiling(averted_outcomes_diff_pax * 0.8 * 0.895),
         total_averted_outcomes = averted_outcomes_cascade_pax + averted_outcomes_diff_vax,
         perc_averted = round(total_averted_outcomes/ total_death * 100, 1),
         strategy = factor(strategy, levels = c("strat3", "strat4", "strat5", "strat6")))
         

