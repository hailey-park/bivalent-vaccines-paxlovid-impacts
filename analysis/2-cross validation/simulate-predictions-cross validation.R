###################################################################################################
#Title: Prediction Simulation
#Author: Hailey Park
#Date: January 10, 2023
###################################################################################################

rm(list=ls())

setwd("/mnt/projects/covid_partners/ucsf_lo/Bivalent Booster Paxlovid Impact Analysis")

#Load in libraries
library(readr)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(tidyr)
library(tibble)
library(reshape2)
library(lubridate)


#Read in dataset
covid_case_data <- read.csv("outcome data clean/covid_data_clean_cases_new.csv")[,-1]
covid_hosp_data <- read.csv("outcome data clean/covid_data_clean_hosp_new.csv")[,-1]
covid_death_data <- read.csv("outcome data clean/covid_data_clean_death_new.csv")[,-1]

#Make dataframes for each outcome (cases, hosp, deaths)
covid_cases <- covid_case_data %>%
  mutate(weeks_since_july2022 = as.integer(floor(difftime(strptime(dtepisode, format="%Y-%m-%d"), 
                                                            strptime("2022-07-23", format="%Y-%m-%d"), 
                                                            units = 'weeks')))) %>%
  group_by(weeks_since_july2022, age_group, boost_2_vax_status_match, .drop = FALSE) %>%
  summarise(num_cases = n()) 

covid_hosp <- covid_hosp_data  %>%
  mutate(weeks_since_july2022 = as.integer(floor(difftime(strptime(dt_hosp, format="%Y-%m-%d"), 
                        strptime("2022-07-23", format="%Y-%m-%d"), 
                        units = 'weeks')))) %>%
  group_by(weeks_since_july2022, age_group, boost_2_vax_status_match, .drop = FALSE) %>%
  summarise(num_hosp = n()) %>% mutate(num_hosp = (num_hosp * (1/0.75)))

covid_death <- covid_death_data %>%
  mutate(weeks_since_july2022 = as.integer(floor(difftime(strptime(death_date_both_or, format="%Y-%m-%d"), 
                                                            strptime("2022-07-23", format="%Y-%m-%d"), 
                                                            units = 'weeks')))) %>%
  group_by(weeks_since_july2022, age_group, boost_2_vax_status_match, .drop = FALSE) %>%
  summarise(num_death = n()) 

# fix n as number of simulations
n <- 1000


#Generate random probabilities
set.seed(42)
p <- matrix(runif(n*3, 0, 1), nrow = n, ncol=3)


#Create regression model to generate predictions
model1 <- glm(num_cases ~ age_group + boost_2_vax_status_match, 
              data = covid_cases %>% filter(weeks_since_july2022 %in% c(0:26)), "quasipoisson")
calibration_period <- covid_cases %>% filter(weeks_since_july2022 %in% c(0:26))
pred <- predict(model1, newdata = calibration_period, se.fit = T)
calibration_period$prediction <- exp(pred$fit)
calibration_period$sd <- exp(pred$se.fit)

sum(calibration_period$prediction)

list <- lapply(p[,1], function(x) qnorm(x, calibration_period$prediction, calibration_period$sd))
  
vax_cases <- do.call(rbind, list)
vax_cases[vax_cases < 0] <- 0

saveRDS(vax_cases, "data/simulated-cases-crossval.RDS")
saveRDS(vax_hosp, "data/simulated-hosp-crossval.RDS")
saveRDS(vax_death, "data/simulated-death-crossval.RDS")
saveRDS(calibration_period %>% select(weeks_since_july2022, age_group, boost_2_vax_status_match), 'data/base-dataframe-crossval.RDS')








