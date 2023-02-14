###################################################################################################
#Title: COVID Vaccine Data Cleaning
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

#Read in csv
covid_vaccine_data <- read.csv('covid19vaccinesadministeredbydemographics.csv')

#Data cleaning
covid_vaccine_data$administered_date <- as.Date(covid_vaccine_data$administered_date)

#Creating bivalent doses distributed for each vaccine strategy
all_ca <- covid_vaccine_data %>%
  filter(demographic_category == 'Age Group',
         administered_date == as.Date('2023-01-23'),
         demographic_value != "Unknown Agegroup")


over_50 <- covid_vaccine_data %>%
  filter(demographic_category == 'Age Group',
         administered_date == as.Date('2023-01-23'),
         demographic_value %in% c("65+", "50-64"))
over_65 <- covid_vaccine_data %>%
  filter(demographic_category == 'Age Group',
         administered_date == as.Date('2023-01-23'),
         demographic_value %in% c("65+"))

prop_over_75 <- 6.1/(6.1 + 9.1) #proportion of 65+ that are 75+
prop_over_50 <- .336 #proportion of CA pop that are 50+
prop_over_65 <- .156 #proportion of CA pop that are 65+ (from KFF website)
unvaccinated_count <- 39240000 - sum(all_ca$cumulative_at_least_one_dose)
unboosted_count <- sum(all_ca$cumulative_fully_vaccinated) - sum(all_ca$cumulative_booster_recip_count)
previously_vacc_count <- sum(all_ca$cumulative_at_least_one_dose)
boosted_count <- sum(all_ca$cumulative_booster_recip_count)


bivalent_booster_doses <- data.frame(strategy = c("strat1", "strat2", "strat3", "strat4", "strat5", "strat6", "strat7", "strat8", "strat9", "strat10", "strat11", "strat12", "strat13", "strat14", "strat15", "strat16"),
                                     bivalent_doses = c(39240000, #strat1
                                                        previously_vacc_count, #strat2
                                                        unvaccinated_count, #strat3
                                                        unboosted_count, #strat4
                                                        sum(over_65$cumulative_at_least_one_dose) * prop_over_75, #strat5
                                                        sum(over_65$cumulative_at_least_one_dose), #strat6
                                                        sum(over_50$cumulative_at_least_one_dose), #strat7
                                                        (39240000 * prop_over_50) - sum(over_50$cumulative_fully_vaccinated), #strat8
                                                        (39240000 * prop_over_65) - sum(over_65$cumulative_fully_vaccinated), #strat9
                                                        ((39240000 * prop_over_65) - (sum(over_65$cumulative_fully_vaccinated))) * prop_over_75, #strat10
                                                        sum(over_50$cumulative_fully_vaccinated) - sum(over_50$cumulative_booster_recip_count), #strat11
                                                        sum(over_65$cumulative_fully_vaccinated) - sum(over_65$cumulative_booster_recip_count), #strat12
                                                        (sum(over_65$cumulative_fully_vaccinated) - sum(over_65$cumulative_booster_recip_count)) * prop_over_75, #strat13
                                                        sum(over_50$cumulative_booster_recip_count), #strat14
                                                        sum(over_65$cumulative_booster_recip_count), #strat15
                                                        sum(over_65$cumulative_booster_recip_count) * prop_over_75)) #strat16


bivalent_booster_clean <- melt(bivalent_booster_doses) %>% rename(total_doses = value) %>% mutate(total_doses = ceiling(total_doses))

write.csv(bivalent_booster_clean , 'bivalent_booster_doses_clean_new.csv')

