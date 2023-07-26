###################################################################################################
#Title: COVID Case Data Cleaning (NEW)
#Author: Hailey Park
#Date: July 10th, 2023
###################################################################################################

rm(list=ls())

setwd("/mnt/projects/covid_partners/ucsf_lo")


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
library(gridExtra)
library(extrafont)
font_import()
loadfonts(device="win") 
fonts()                   



#Loading in datasets
covid_data1 <- fread("covid_registry_2023-05-09.csv", 
                     select = c("age", 
                                "dtepisode",
                                "raceeth",
                                "sex2",
                                "cste_def",
                                "hospitalized",
                                "died",
                                "test_type",
                                "boost_2_vax_status_match"
                     ), 
                     showProgress = FALSE)


column_names <- as.data.frame(colnames(read.csv(file = "covid_registry_2023-05-09.csv", nrows = 100)))

####################################################################################################
#COVID Data Cleaning

covid_data_cleaning <- covid_data1

#Converting dated variables to Date objects
covid_data_cleaning$dtepisode <- as.Date(covid_data_cleaning$dtepisode)
covid_data_cleaning$dtdeath <- as.Date(covid_data_cleaning$dtdeath)

#Only keeping confirmed cases and reasonable ages (0-110 yrs old) and data from the last 6 mo (7/23/22 - 1/23/23)
covid_data_cleaning <- covid_data_cleaning %>%
  filter(age >= 0,
         age <= 110,
         cste_def == "Confirmed",
         (dtepisode >= '2022-07-23' & dtepisode <= '2023-01-23'))

#Creating hospitalizations date variable
covid_data_cleaning <- covid_data_cleaning %>%
  mutate(dt_hosp = if_else(hospitalized == "Y",
                           dtepisode,
                           NULL))


#Creating death date variables
covid_data_cleaning <- covid_data_cleaning %>%
  mutate(dt_death = if_else(died == 'Y',
                            dtepisode,
                            NULL))


#Creating an age-group variable
covid_data_cleaning <- covid_data_cleaning %>% 
  mutate(
    # Create categories
    age_group = dplyr::case_when(
      age <= 17            ~ "0-17 years",
      age > 17 & age <= 49 ~ "18-49 years",
      age > 49 & age <= 64 ~ "50-64 years",
      age > 64 & age <= 74 ~ "65-74 years",
      age > 74 & age <= 84 ~ "75-84 years",
      age > 84 ~ "85+ years"
    ),
    # Convert to factor
    age_group = factor(
      age_group,
      level = c("0-17 years", "18-49 years","50-64 years","65-74 years" ,"75-84 years", "85+ years")
    )
  )

#Relabeling the vaccination status values
covid_data_cleaning <- covid_data_cleaning %>%
  mutate(boost_2_vax_status_match = recode(boost_2_vax_status_match, 
                                           'partial_vax' = 'Primary Series',
                                           'full_vax' = 'Primary Series',
                                           'boosted - one dose' = 'Boosted (1 dose)',
                                           'boosted - two dose' = 'Boosted (2 doses)'),
         boost_2_vax_status_match = replace(boost_2_vax_status_match, 
                                            boost_2_vax_status_match == "",
                                            'Unvaccinated')) 

#Re-ordering th vaccination status values
covid_data_cleaning$boost_2_vax_status_match <- factor(covid_data_cleaning$boost_2_vax_status_match, 
                                                       levels = c("Unvaccinated", 
                                                                  "Primary Series", 
                                                                  "Boosted (1 dose)",
                                                                  "Boosted (2 doses)"))
#Save as .csv file
write.csv(covid_data_cleaning %>% filter(!is.na(dt_death)), 'covid_data_clean_death_final.csv')
write.csv(covid_data_cleaning %>% filter(!is.na(dt_hosp)), 'covid_data_clean_hosp_final.csv')
write.csv(covid_data_cleaning, 'covid_data_clean_cases_final.csv')
####################################################################################################

inspection <- covid_data_cleaning %>% #filter(!is.na(death_date_both_or)) %>%
  group_by(sex2) %>% summarise(total_outcomes = n()) %>% 
  mutate(perc_outcomes = round(total_outcomes/sum(total_outcomes) * 100, 2))


sum_outcomes <- sum(inspection$total_outcomes)
