###################################################################################################
#Title: Nirmatrelvir-Ritonavir Regression Analysis
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
covid_case_data <- read.csv("outcome data clean/covid_data_clean_cases_new.csv")[,-1]
covid_hosp_data <- read.csv("outcome data clean/covid_data_clean_hosp_new.csv")[,-1]
covid_death_data <- read.csv("outcome data clean/covid_data_clean_death_new.csv")[,-1]

#Make dataframes for each outcome (cases, hosp, deaths)
covid_cases <- covid_case_data %>%
  group_by(dtepisode, age_group, boost_2_vax_status_match, .drop = FALSE) %>%
  summarise(num_cases = n())


covid_hosp <-covid_hosp_data  %>%
  group_by(dt_hosp, age_group, boost_2_vax_status_match, .drop = FALSE) %>%
  summarise(num_hosp = n()) %>% mutate(num_hosp = (num_hosp * (1/0.75))) 


covid_death <- covid_death_data %>%
  group_by(death_date_both_or, age_group, boost_2_vax_status_match, .drop = FALSE) %>%
  summarise(num_death = n())

pax_death <- merge(covid_cases, covid_death, by.x = c("dtepisode", "age_group", "boost_2_vax_status_match"), 
                   by.y = c("death_date_both_or", "age_group", "boost_2_vax_status_match"), all.x = TRUE) %>% replace(is.na(.), 0)

pax_hosp <- merge(covid_cases, covid_hosp, by.x = c("dtepisode", "age_group", "boost_2_vax_status_match"), 
                  by.y = c("dt_hosp", "age_group", "boost_2_vax_status_match"), all.x = TRUE)


pe_estimates <- data.frame(baseline_vacc = c("Unvaccinated", "Vaccinated"),
                           pe_hosps = c(.89, .40),
                           pe_deaths = c(.89, .71))


#turning Date object 'dtepisode' variables to 'weeks_since_july2022'
pax_hosp$weeks_since_july2022 <- as.integer(floor(difftime(strptime(pax_hosp$dtepisode, format="%Y-%m-%d"),
                                                           strptime("2022-07-23", format="%Y-%m-%d"),
                                                           units = 'weeks')))

pax_death$weeks_since_july2022 <- as.integer(floor(difftime(strptime(pax_death$dtepisode, format="%Y-%m-%d"),
                                                            strptime("2022-07-23", format="%Y-%m-%d"),
                                                            units = 'weeks')))

###################################################################################################
#Regression Model: Simple Model using only age_group and vaccination status as predictors

#Model Calibration (on 6 months of data)
model1 <- glm(num_hosp ~ age_group + boost_2_vax_status_match, 
              data = pax_hosp %>% filter(weeks_since_july2022 %in% c(0:26)), "quasipoisson")
summary(model1)

#Model Calibration Check (Predicting on 6 months of data that it was calibrated on)
calibration_period <- pax_hosp %>% filter(weeks_since_july2022 %in% c(0:26))
pred <- predict(model1, newdata = calibration_period, se.fit = T)
calibration_period$prediction_hosp <- exp(pred$fit)

sum(calibration_period$num_hosp)

###################################################################################################
#Calculating Averted Outcomes

#Creating predictions df
predictions_pax_df <- merge(calibration_period %>% mutate(vacc_status = ifelse(boost_2_vax_status_match %in% c("Primary Series", "Boosted (1 dose)", "Boosted (2 doses)"),
                                                                               "Vaccinated",
                                                                               "Unvaccinated")) %>% group_by(age_group, boost_2_vax_status_match, vacc_status) %>% 
                              summarise(num_cases = sum(num_cases), 
                                        predicted_cases = sum(prediction_cases), 
                                        num_hosp = sum(num_hosp),
                                        predicted_hosp = sum(prediction_hosp)),
                            pe_estimates, by.x = "vacc_status", by.y = "baseline_vacc", all.x = TRUE)



#######################################
#Primary Strategies

#Strategy 1: Giving Paxlovid to those 18+ years (everyone)
strat1 <- predictions_pax_df %>% mutate(averted_hosp = case_when((age_group %in% c("18-49 years", "50-64 years", "65-74 years", "75-84 years", "85+ years"))  ~ predicted_hosp * pe_hosps,
                                                                 TRUE ~ 0),
                                        target_group_cases = case_when((age_group %in% c("18-49 years", "50-64 years", "65-74 years", "75-84 years", "85+ years")) ~ predicted_cases,
                                                                       TRUE ~ 0))
strat1$strategy <- rep("strat1", nrow(strat1))

#Strategy 2: Giving Paxlovid to those 50+ years (everyone)
strat2 <- predictions_pax_df %>% mutate(averted_hosp = case_when((age_group %in% c("50-64 years", "65-74 years", "75-84 years", "85+ years"))  ~ predicted_hosp * pe_hosps,
                                                                 TRUE ~ 0),
                                        target_group_cases = case_when((age_group %in% c("50-64 years", "65-74 years", "75-84 years", "85+ years")) ~ predicted_cases,
                                                                       TRUE ~ 0))
strat2$strategy <- rep("strat2", nrow(strat2))

#Strategy 3: Giving Paxlovid to those 50+ years (unvaccinated and comorbidities)
strat3 <- predictions_pax_df %>% mutate(hosp_unvaccinated = ifelse(age_group == "50-64 years" & vacc_status %in% c("Unvaccinated"), predicted_hosp, 0),
                                        hosp_comorbidities = ifelse(age_group == "50-64 years", predicted_hosp * 0.33, 0),
                                        hosp_both = ifelse(age_group == "50-64 years" & vacc_status %in% c("Unvaccinated"), predicted_hosp * 0.33, 0),
                                        hosp_either = hosp_unvaccinated + hosp_comorbidities - hosp_both,
                                        cases_unvaccinated = ifelse(age_group == "50-64 years" & vacc_status %in% c("Unvaccinated"), predicted_cases, 0),
                                        cases_comorbidities = ifelse(age_group == "50-64 years", predicted_cases * 0.33, 0),
                                        cases_both = ifelse(age_group == "50-64 years" & vacc_status %in% c("Unvaccinated"), predicted_cases * 0.33, 0),
                                        cases_either = cases_unvaccinated + cases_comorbidities - cases_both,
                                        averted_hosp = case_when((age_group == "50-64 years") ~ hosp_either * pe_hosps,
                                                                 (age_group %in% c("65-74 years", "75-84 years", "85+ years"))  ~ predicted_hosp * pe_hosps,
                                                                 TRUE ~ 0),
                                        target_group_cases = case_when((age_group == "50-64 years") ~ cases_either,
                                                                       (age_group %in% c("65-74 years", "75-84 years", "85+ years")) ~ predicted_cases,
                                                                       TRUE ~ 0))

strat3$strategy <- rep("strat3", nrow(strat3))


#Strategy 4: Giving Paxlovid to the those 65+ years
strat4 <- predictions_pax_df %>% mutate(averted_hosp = ifelse((age_group %in% c("65-74 years", "75-84 years", "85+ years")),
                                                              predicted_hosp * pe_hosps, 0),
                                        target_group_cases = ifelse((age_group %in% c("65-74 years", "75-84 years", "85+ years")),
                                                                    predicted_cases, 0))

strat4$strategy <- rep("strat4", nrow(strat4))

#Strategy 5: Giving Paxlovid to the those 75+ years
strat5 <- predictions_pax_df %>% mutate(averted_hosp = ifelse((age_group %in% c("75-84 years", "85+ years")),
                                                              predicted_hosp * pe_hosps, 0),
                                        target_group_cases = ifelse((age_group %in% c("75-84 years", "85+ years")),
                                                                    predicted_cases, 0) )

strat5$strategy <- rep("strat5", nrow(strat5))  


combined_strat <- bind_rows(strat1, strat2, strat3, strat4, strat5)

#######################################
#Additional Paxlovid Strategies

#Strategy 1: Giving Paxlovid to those Unvaccinated (65+ years)
strat1 <- predictions_pax_df %>% mutate(averted_hosp = ifelse((age_group %in% c("65-74 years", "75-84 years", "85+ years") & boost_2_vax_status_match %in% c("Unvaccinated")), predicted_hosp * pe_hosps, 0),
                                        target_group_cases = ifelse((age_group %in% c("65-74 years", "75-84 years", "85+ years") & boost_2_vax_status_match %in% c("Unvaccinated")), predicted_cases, 0))

strat1$strategy <- rep("strat1", nrow(strat1))

#Strategy 2: Giving Paxlovid to those Unvaccinated (75+ years)
strat2 <- predictions_pax_df %>% mutate(averted_hosp = ifelse((age_group %in% c("75-84 years", "85+ years") & boost_2_vax_status_match %in% c("Unvaccinated")), predicted_hosp * pe_hosps, 0),
                                        target_group_cases = ifelse((age_group %in% c("75-84 years", "85+ years") & boost_2_vax_status_match %in% c("Unvaccinated")), predicted_cases, 0))

strat2$strategy <- rep("strat2", nrow(strat2))


#Strategy 3:Giving Paxlovid to those Primary Series only (65+ years)
strat3 <- predictions_pax_df %>% mutate(averted_hosp = ifelse((age_group %in% c("65-74 years", "75-84 years", "85+ years") & boost_2_vax_status_match %in% c("Primary Series")), predicted_hosp * pe_hosps, 0),
                                        target_group_cases = ifelse((age_group %in% c("65-74 years", "75-84 years", "85+ years") & boost_2_vax_status_match %in% c("Primary Series")), predicted_cases, 0))

strat3$strategy <- rep("strat3", nrow(strat3))

#Strategy 4: Giving Paxlovid to those Primary Series only (75+ years)
strat4 <- predictions_pax_df %>% mutate(averted_hosp = ifelse((age_group %in% c("75-84 years", "85+ years") & boost_2_vax_status_match %in% c("Primary Series")), predicted_hosp * pe_hosps, 0),
                                        target_group_cases = ifelse((age_group %in% c("75-84 years", "85+ years") & boost_2_vax_status_match %in% c("Primary Series")), predicted_cases, 0))

strat4$strategy <- rep("strat4", nrow(strat4))

#Strategy 5: Giving Paxlovid to those Boosted only (65+ years)
strat5 <- predictions_pax_df %>% mutate(averted_hosp = ifelse((age_group %in% c("65-74 years", "75-84 years", "85+ years") & boost_2_vax_status_match %in% c("Boosted (1 dose)", "Boosted (2 doses)")), predicted_hosp * pe_hosps, 0),
                                        target_group_cases = ifelse((age_group %in% c("65-74 years", "75-84 years", "85+ years") & boost_2_vax_status_match %in% c("Boosted (1 dose)", "Boosted (2 doses)")), predicted_cases, 0))

strat5$strategy <- rep("strat5", nrow(strat5))

#Strategy 6: Giving Paxlovid to those Boosted only (75+ years)
strat6 <- predictions_pax_df %>% mutate(averted_hosp = ifelse((age_group %in% c("75-84 years", "85+ years") & boost_2_vax_status_match %in% c("Boosted (1 dose)", "Boosted (2 doses)")), predicted_hosp * pe_hosps, 0),
                                        target_group_cases = ifelse((age_group %in% c("75-84 years", "85+ years") & boost_2_vax_status_match %in% c("Boosted (1 dose)", "Boosted (2 doses)")), predicted_cases, 0))

strat6$strategy <- rep("strat6", nrow(strat6))

combined_strat <- bind_rows(strat1, strat2, strat3, strat4, strat5, strat6)


#######################################
combined_strat_outcomes <- combined_strat %>% group_by(strategy) %>% summarise(averted_outcomes_max = ceiling(sum(averted_hosp)),
                                                                               total_hosp = ceiling(sum(predicted_hosp)),
                                                                               target_group_cases = ceiling(sum(target_group_cases))) %>%
  mutate(averted_outcomes_cascade = ceiling(averted_outcomes_max * 0.8 * 0.895),
         treated_cascade = ceiling(target_group_cases * 0.8 *0.895),
         current_averted = ceiling(averted_outcomes_cascade * 0.3),
         current_treated = ceiling(treated_cascade * 0.3),
         actual_averted = averted_outcomes_cascade - current_averted,
         actual_treated = treated_cascade - current_treated,
         NNT = ceiling(actual_treated/actual_averted),
         perc_averted = round(actual_averted / total_hosp * 100, 1),
         strategy = factor(strategy, levels = c("strat1", "strat4", "strat2", "strat5", "strat3", "strat6")))



