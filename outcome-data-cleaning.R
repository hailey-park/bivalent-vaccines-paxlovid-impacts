###################################################################################################
#Title: COVID Case Data Cleaning
#Author: Hailey Park
#Date: January 30th, 2023
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
covid_data1 <- fread("covid_registry_2023-01-30.csv", 
                     select = c("age", 
                                "cntyfips",
                                "county", 
                                "ctract",
                                "dtdeath",
                                "dtepisode",
                                "raceeth",
                                "sex2",
                                "cste_def",
                                "hospitalized",
                                "died",
                                "test_type",
                                "boost_2_vax_status_match",
                                "postvax_case_match",
                                "ncovpuihospdtladmitdt_1",
                                "ncovpuihospdtladmitdt_2",
                                "ncovpuihospdtladmitdt_3",
                                "dtlabresult"
                     ), 
                     showProgress = FALSE)


column_names <- as.data.frame(colnames(read.csv(file = "covid_registry_2023-01-30.csv", nrows = 100)))

####################################################################################################
#COVID Data Cleaning

covid_data_cleaning <- covid_data1

#Converting dated variables to Date objects
covid_data_cleaning$dtepisode <- as.Date(covid_data_cleaning$dtepisode)
covid_data_cleaning$dtdeath <- as.Date(covid_data_cleaning$dtdeath)
covid_data_cleaning$ncovpuihospdtladmitdt_1 <- as.Date(covid_data_cleaning$ncovpuihospdtladmitdt_1)
covid_data_cleaning$ncovpuihospdtladmitdt_2 <- as.Date(covid_data_cleaning$ncovpuihospdtladmitdt_2)
covid_data_cleaning$ncovpuihospdtladmitdt_3 <- as.Date(covid_data_cleaning$ncovpuihospdtladmitdt_3)

#Only keeping confirmed cases and reasonable ages (0-110 yrs old) and data from the last 6 mo (3/1/22 - 9/8/22)
covid_data_cleaning <- covid_data_cleaning %>%
  filter(age >= 0,
         age <= 110,
         cste_def == "Confirmed",
         (dtepisode >= '2022-01-23' & dtepisode <= '2023-01-23'))

#Creating hospitalizations variables to double check counts
covid_data_cleaning <- covid_data_cleaning %>%
  mutate(hospitalized_using_dates = ifelse((!is.na(ncovpuihospdtladmitdt_1) | !is.na(ncovpuihospdtladmitdt_2) | !is.na(ncovpuihospdtladmitdt_3)),
                                           "Y",
                                           "N"),
         both_hospitalized_checks = ifelse(hospitalized == "Y" & hospitalized_using_dates == "Y",
                                           "Y",
                                           "N")) 
covid_data_cleaning$earliest_hospitalized_dt <- do.call(pmin, c(covid_data_cleaning %>%
                                                                  select(ncovpuihospdtladmitdt_1, ncovpuihospdtladmitdt_2, ncovpuihospdtladmitdt_3),
                                                                na.rm = T))

covid_data_cleaning <- covid_data_cleaning %>%
  mutate(dt_hosp = if_else(hospitalized == "Y" & is.na(earliest_hospitalized_dt),
                           dtepisode,
                           if_else(!is.na(earliest_hospitalized_dt),
                                   dtepisode, #####vhange back to earliest_hosp
                                   NULL)))


#Creating death variables to double check counts
covid_data_cleaning <- covid_data_cleaning %>%
  mutate(death_date_both_and = if_else(died == 'Y',
                                       dtdeath,
                                       NULL),
         death_date_both_or = if_else(died == 'Y' & is.na(dtdeath),
                                      dtepisode,
                                      if_else(!is.na(dtdeath),
                                              dtepisode, #######change back to dtdeath
                                              NULL)))


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

inspection <- covid_data_cleaning %>% group_by (boost_2_vax_status_match) %>%
  summarise(total = n())

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
inspection <- covid_data_cleaning %>% filter(!is.na(death_date_both_or))
#Save as .csv file
write.csv(covid_data_cleaning %>% filter(!is.na(death_date_both_or)), 'covid_data_clean_death_crossval-6mo.csv')
write.csv(covid_data_cleaning %>% filter(!is.na(dt_hosp)), 'covid_data_clean_hosp_crossval-6mo.csv')
write.csv(covid_data_cleaning, 'covid_data_clean_cases_crossval-6mo.csv')
####################################################################################################

inspection <- covid_data_cleaning %>% filter(!is.na(dt_hosp),
                                             age_group %in% c("65-74 years", "75-84 years", "85+ years")) %>%
  group_by(boost_2_vax_status_match) %>% summarise(total_outcomes = n()) %>% 
  mutate(perc_outcomes = round(total_outcomes/sum_outcomes * 100, 2))


sum_outcomes <- sum(inspection$total_outcomes)
