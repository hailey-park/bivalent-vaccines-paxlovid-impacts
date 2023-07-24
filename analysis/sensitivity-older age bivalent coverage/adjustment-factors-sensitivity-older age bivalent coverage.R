###################################################################################################
#Title: Adjustment Factors for Bivalent Coverage
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


#Load in dataset
vaccination_data <- read.csv("covid19vaccinesadministeredbydemographics.csv")
vaccination_data$administered_date <- as.Date(vaccination_data$administered_date)
vaccination_data$demographic_value <- as.character(vaccination_data$demographic_value)

#Calculate cumulative bivalent coverage and prevalence of vaccine groups over time
bivalent_coverage <- vaccination_data %>%
  filter(demographic_category == 'Age Group',
         demographic_value != "Unknown Agegroup",
         administered_date > as.Date("2022-06-30"),
         administered_date < as.Date("2023-02-01")) %>%
  group_by(month = lubridate::floor_date(administered_date, "month"), demographic_value) %>%
  summarise(biv_booster_eligible = max(bivalent_booster_eligible_population),
            total_bivalent_cumulative = max(cumulative_bivalent_booster_recip_count),
            total_booster_cumulative = max(cumulative_booster_recip_count) - total_bivalent_cumulative,
            fully_vaccinated = max(cumulative_fully_vaccinated) - total_booster_cumulative) %>%
  mutate(age_group = ifelse(demographic_value %in% c("Under 5", "5-11", "12-17"), "0-17", demographic_value)) %>% 
  group_by(month, age_group) %>% summarise(bivalent = sum(total_bivalent_cumulative),
                                           bivalent_eligible_pop = sum(biv_booster_eligible),
                                           booster = sum(total_booster_cumulative),
                                           fully_vaccinated = sum(fully_vaccinated)) %>%
  mutate(percent_bivalent_cum = round(bivalent/bivalent_eligible_pop, 2))

#Create pop df and merge with bivalent coverage df to get prevalence of each vaccine group
pop_by_age <- data.frame(age_group = c("0-17", "18-49", "50-64", "65+"),
                         total_pop = c(ceiling(39240000 * .224),      #got these age-bin multipliers from Census/Statista
                                       ceiling(39240000 * .446),
                                       ceiling(39240000 * .18),
                                       ceiling(39240000 * .15)))

bivalent_coverage_prev <- merge(bivalent_coverage, pop_by_age, by.x = "age_group", by.y = "age_group", all.x = TRUE) %>%
  mutate(prev_boost = round(booster/total_pop, 2),
         prev_prim = round(fully_vaccinated/total_pop, 2),
         prev_boost2 = case_when(age_group %in% c("50-64", "65+") ~ prev_boost * 0.4,
                                 age_group == "18-49" ~ prev_boost * 0.3,
                                 age_group == "0-17" ~ prev_boost * 0.1,
                                 TRUE ~ 0),
         prev_boost1 = prev_boost - prev_boost2,
         prev_unvax = 1 - prev_prim - prev_boost1 - prev_boost2,
         percent_bivalent_cum_v2 = round(bivalent/total_pop, 2)) %>%
  dplyr::select(age_group, month, percent_bivalent_cum, percent_bivalent_cum_v2, prev_unvax, prev_prim, prev_boost1, prev_boost2)

#Add the oldest age groups (65-74, 75-84, 85+)
oldest_age <- do.call("rbind", replicate(3, bivalent_coverage_prev %>% 
                                           filter(age_group == "65+"), simplify = FALSE))

oldest_age$age_group <- c(rep("65-74", 7), rep("75-84", 7), rep("85+", 7))

oldest_age_adjusted <- oldest_age %>%
  mutate(percent_bivalent_cum_v2 = case_when(age_group == "85+" ~ percent_bivalent_cum_v2 + 0.1,
                                              age_group == "75-84" ~ percent_bivalent_cum_v2 + 0.05,
                                              age_group == "65-74" ~ ((percent_bivalent_cum_v2 - (percent_bivalent_cum_v2 + 0.05) * (4.4/15.2) - (percent_bivalent_cum_v2 + 0.1) * (1.7/15.2)) * (15.2/9.1))))

bivalent_coverage_with_old <- rbind(bivalent_coverage_prev, oldest_age_adjusted) %>% filter(age_group != "65+") %>%
  mutate(total_prev = prev_unvax + prev_prim + prev_boost1 + prev_boost2)


#Calculate adjustment factors
adjustment_factor_matrix <- bivalent_coverage_with_old %>% 
  mutate(perc_prim = case_when(age_group %in% c("50-64", "65-74", "75-84", "85+") ~ 0.2,
                               age_group == "18-49" ~ 0.09,
                               age_group == "0-17" ~ 0.14,
                               TRUE ~ 0),
         perc_boost1 = case_when(age_group %in% c("50-64", "65-74", "75-84", "85+") ~ 0.4,
                                 age_group == "18-49" ~ 0.27,
                                 age_group == "0-17" ~ 0.28,
                                 TRUE ~ 0),
         perc_boost2 = case_when(age_group %in% c("50-64", "65-74", "75-84", "85+") ~ 0.52,
                                 age_group == "18-49" ~ 0,
                                 age_group == "0-17" ~ 0,
                                 TRUE ~ 0),
         prev_biv_prim = perc_prim * prev_prim,
         prev_biv_boost1 = perc_boost1 * prev_boost1,
         prev_biv_boost2 = perc_boost2 * prev_boost2,
         solve_x = percent_bivalent_cum_v2/(prev_biv_prim + prev_biv_boost1 + prev_biv_boost2),
         prev_biv_prim_abs = prev_biv_prim * solve_x,
         prev_biv_boost1_abs = prev_biv_boost1 * solve_x,
         prev_biv_boost2_abs = prev_biv_boost2 * solve_x,
         prev_biv_prim_rel = prev_biv_prim_abs / prev_prim,
         prev_biv_boost1_rel = prev_biv_boost1_abs / prev_boost1,
         prev_biv_boost2_rel = prev_biv_boost2_abs / prev_boost2,
         ve_prim = 0.41, #CHANGE HERE
         ve_boost1 = 0.26, #CHANGE HERE
         ve_boost2 = 0.40, #CHANGE HERE
         adj_prim = round(1/(1 - (prev_biv_prim_rel * ve_prim)), 2),
         adj_boost1 = round(1/(1 - (prev_biv_boost1_rel * ve_boost1)), 2),
         adj_boost2 = round(1/(1 - (prev_biv_boost2_rel * ve_boost2)), 2)) %>%
  arrange(month)

#cleaning up df
adjustment_factor_matrix_melt <- melt(adjustment_factor_matrix %>%
                                        select(age_group, month, adj_prim, adj_boost1, adj_boost2), id = c("age_group", "month")) %>%
  mutate(age_group = paste0(age_group, " years"),
         adj_factor = value,
         boost_2_vax_status_match = case_when(variable == "adj_prim" ~ "Primary Series",
                                              variable == "adj_boost1" ~ "Boosted (1 dose)",
                                              variable == "adj_boost2" ~ "Boosted (2 doses)",
                                              TRUE ~ "NA"))
adjustment_factor_clean <- rbind(adjustment_factor_matrix_melt,
                                 (adjustment_factor_matrix_melt %>% filter(variable == "adj_prim")%>%
                                    mutate(adj_factor = 1,
                                           boost_2_vax_status_match = "Unvaccinated"))) %>% select(-c(variable, value))

write.csv(adjustment_factor_clean, "outcome data clean/adjustment_factor_matrix_cases-sensitivity-oldest age bivalent coverage.csv")

###########################################################
#Bivalent coverage by age group and vaccine status

bivalent_coverage_age_vax <- melt(adjustment_factor_matrix %>%
                                    select(age_group, month, prev_biv_prim_rel, prev_biv_boost1_rel, prev_biv_boost2_rel), id = c("age_group", "month")) %>%
  mutate(age_group = paste0(age_group, " years"),
         bivalent_coverage = value,
         boost_2_vax_status_match = case_when(variable == "prev_biv_prim_rel" ~ "Primary Series",
                                              variable == "prev_biv_boost1_rel" ~ "Boosted (1 dose)",
                                              variable == "prev_biv_boost2_rel" ~ "Boosted (2 doses)",
                                              TRUE ~ "NA"))
bivalent_coverage_age_vax_clean <- rbind(bivalent_coverage_age_vax,
                                         (bivalent_coverage_age_vax %>% filter(variable == "prev_biv_prim_rel")%>%
                                            mutate(bivalent_coverage = 0,
                                                   boost_2_vax_status_match = "Unvaccinated"))) %>% 
  select(-c(variable, value)) %>%
  filter(month == "2022-11-01") %>%
  select(age_group, boost_2_vax_status_match, bivalent_coverage)


write.csv(bivalent_coverage_age_vax_clean, "outcome data clean/bivalent_coverage_age_vax-sensitivity-oldest age bivalent coverage.csv")







