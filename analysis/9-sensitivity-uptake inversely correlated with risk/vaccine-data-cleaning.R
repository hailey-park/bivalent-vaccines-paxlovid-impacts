###################################################################################################
#Title: COVID Vaccine Data Cleaning (incorporating reverse uptake)
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
over_65 <- covid_vaccine_data %>%
  filter(demographic_category == 'Age Group',
         administered_date == as.Date('2023-01-23'),
         demographic_value %in% c("65+"))

under_65 <- covid_vaccine_data %>%
  filter(demographic_category == 'Age Group',
         administered_date == as.Date('2023-01-23'),
         demographic_value != "65+")

age_50_64 <- covid_vaccine_data %>%
  filter(demographic_category == 'Age Group',
         administered_date == as.Date('2023-01-23'),
         demographic_value == "50-64")


prop_over_65 <- .156 #proportion of CA pop that are 65+ (from KFF website)
prop_over_75 <- 6.1/(6.1 + 9.1) #proportion of 65+ that are 75+

total_over_65 <- 39240000 * prop_over_65
unvaccinated_count_over_65 <- total_over_65 - sum(over_65$cumulative_fully_vaccinated)
previously_vacc_count_over_65 <- sum(over_65$cumulative_fully_vaccinated)
bivalent_count_over_65 <- sum(over_65$cumulative_bivalent_booster_recip_count)
primary_series_count_over_65 <- (sum(over_65$cumulative_fully_vaccinated) - (sum(over_65$cumulative_booster_recip_count) - bivalent_count_over_65))
#'primary_series_count' already subtracting bivalent count because general cumulative 'booster' includes bivalent counts

total_under_65 <- 39240000 * (1 - prop_over_65)
unvaccinated_count_under_65 <- total_under_65 - sum(under_65$cumulative_fully_vaccinated)
previously_vacc_count_under_65 <- sum(under_65$cumulative_fully_vaccinated)
bivalent_count_under_65 <- sum(under_65$cumulative_bivalent_booster_recip_count)
primary_series_count_under_65 <- (sum(under_65$cumulative_fully_vaccinated) - (sum(under_65$cumulative_booster_recip_count) - bivalent_count_under_65))

previously_vacc_count_50_64 <- sum(age_50_64$cumulative_fully_vaccinated)
bivalent_count_50_64 <- sum(age_50_64$cumulative_bivalent_booster_recip_count)


#Strategy 1: Everyone
#Strategy 2: Previously Vaccinated (primary series, boosted)
#Strategy 3: Unvaccinated
#Strategy 4: Primary Series Only
#Strategy 5: 75+ years (excluding unvaccinated)
#Strategy 6: 65+ years (excluding unvaccinated)
#Strategy 7: 50+ years (excluding unvaccinated)

strat1 <- (total_over_65*0.5 - bivalent_count_over_65) + (total_under_65 * 0.7 - bivalent_count_under_65)
strat2 <- (previously_vacc_count_over_65 * 0.5 - bivalent_count_over_65) + (previously_vacc_count_under_65  * 0.7 - bivalent_count_under_65)
strat3 <- unvaccinated_count_over_65 * (0.5) + unvaccinated_count_under_65 * 0.7
strat4 <- (primary_series_count_over_65 * 0.5 - bivalent_count_over_65) + (primary_series_count_under_65* 0.7 - bivalent_count_under_65)
strat5 <- (previously_vacc_count_over_65*(0.5) - bivalent_count_over_65) * prop_over_75 
strat6 <- (previously_vacc_count_over_65*(0.5) - bivalent_count_over_65) 
strat7 <- (previously_vacc_count_over_65*(0.5) - bivalent_count_over_65) + (previously_vacc_count_50_64*(0.7) - bivalent_count_50_64)



bivalent_booster_doses <- data.frame(strategy = c("strat1", "strat2", "strat3", "strat4", "strat5", "strat6", "strat7"),
                                     bivalent_doses = c(strat1, strat2, strat3,strat4,strat5, strat6, strat7))

bivalent_booster_clean <- melt(bivalent_booster_doses) %>% rename(total_doses = value) %>% mutate(total_doses = ceiling(total_doses))
  
write.csv(bivalent_booster_clean , 'final-code/outcome data clean/bivalent_booster_doses_clean_sensitivity_reverse_uptake.csv')








