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
previously_vacc_count <- sum(all_ca$cumulative_at_least_one_dose)
bivalent_count <- sum(all_ca$cumulative_bivalent_booster_recip_count)
bivalent_eligible <- sum(all_ca$bivalent_booster_eligible_population)
booster_count <- sum(all_ca$cumulative_booster_recip_count)
primary_series_count <- (sum(all_ca$cumulative_fully_vaccinated) - (booster_count - bivalent_count)) 
#'primary_series_count' already subtracting bivalent count because general cumulative 'booster' includes bivalent counts


#Strategy 1: Everyone
#Strategy 2: Previously Vaccinated (primary series, boosted)
#Strategy 3: Unvaccinated
#Strategy 4: Primary Series Only
#Strategy 5: 75+ years (excluding unvaccinated)
#Strategy 6: 65+ years (excluding unvaccinated)
#Strategy 7: 50+ years (excluding unvaccinated)

strat1 <- (39240000*0.7) - bivalent_count
strat2 <- (previously_vacc_count * 0.7) - bivalent_count
strat3 <- unvaccinated_count * 0.7
strat4 <- (primary_series_count * 0.7 - bivalent_count)
strat5 <- sum(over_65$cumulative_at_least_one_dose * 0.7 - over_65$cumulative_bivalent_booster_recip_count) * prop_over_75
strat6 <- sum(over_65$cumulative_at_least_one_dose * 0.7 - over_65$cumulative_bivalent_booster_recip_count)
strat7 <- sum(over_50$cumulative_at_least_one_dose * 0.7 - over_50$cumulative_bivalent_booster_recip_count)

#Strategy 8: Unvaccinated & 50+ years
#Strategy 9: Unvaccinated & 65+ years
#Strategy 10: Unvaccinated & 75+ years
#Strategy 11: Primary Series Only & 50+ years
#Strategy 12: Primary Series Only & 65+ years
#Strategy 13: Primary Series Only & 75+ years
#Strategy 14: Monovalently Boosted Only & 50+ years
#Strategy 15: Monovalently Boosted Only & 65+ years
#Strategy 16: Monovalently Boosted Only & 75+ years

strat8 <- ((39240000 * prop_over_50) - sum(over_50$cumulative_fully_vaccinated)) * 0.7
strat9 <- ((39240000 * prop_over_65) - sum(over_65$cumulative_fully_vaccinated)) * 0.7
strat10 <- ((39240000 * prop_over_65) - (sum(over_65$cumulative_fully_vaccinated))) * prop_over_75 * 0.7
strat11 <- (sum(over_50$cumulative_fully_vaccinated) - (sum(over_50$cumulative_booster_recip_count) - sum(over_50$cumulative_bivalent_booster_recip_count))) * 0.7  - sum(over_50$cumulative_bivalent_booster_recip_count)#cumulative_booster_recip_count also includes bivalent counts
strat12 <- (sum(over_65$cumulative_fully_vaccinated) - (sum(over_65$cumulative_booster_recip_count) - sum(over_65$cumulative_bivalent_booster_recip_count))) * 0.7  - sum(over_65$cumulative_bivalent_booster_recip_count)
strat13 <- ((sum(over_65$cumulative_fully_vaccinated) - (sum(over_65$cumulative_booster_recip_count) - sum(over_65$cumulative_bivalent_booster_recip_count))) * 0.7  - sum(over_65$cumulative_bivalent_booster_recip_count)) * prop_over_75
strat14 <- sum(over_50$cumulative_booster_recip_count)* 0.7 - sum(over_50$cumulative_bivalent_booster_recip_count)
strat15 <- sum(over_65$cumulative_booster_recip_count)* 0.7 - sum(over_65$cumulative_bivalent_booster_recip_count)
strat16 <- (sum(over_65$cumulative_booster_recip_count) * 0.7 - sum(over_65$cumulative_bivalent_booster_recip_count)) * prop_over_75



bivalent_booster_doses <- data.frame(strategy = c("strat1", "strat2", "strat3", "strat4", "strat5", "strat6", "strat7", "strat8", "strat9", "strat10", "strat11", "strat12", "strat13", "strat14", "strat15", "strat16"),
                                     bivalent_doses = c(strat1, strat2, strat3,strat4,strat5, strat6, strat7,strat8, strat9, strat10, strat11, strat12, strat13, strat14, strat15, strat16))

bivalent_booster_clean <- melt(bivalent_booster_doses) %>% rename(total_doses = value) %>% mutate(total_doses = ceiling(total_doses))
  
write.csv(bivalent_booster_clean , 'outcome data clean/bivalent_booster_doses_clean_70maxuptake_v2.csv')
