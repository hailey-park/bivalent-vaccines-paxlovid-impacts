###################################################################################################
#Title: Appendix Figures
#Author: Hailey Park
#Date: July 6, 2023
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
library(extrafont)
           

#########################################################
#Case Fatality and Case Hospitalization Rates Figure

#Read in dataset
covid_case_data <- read.csv("outcome data clean/covid_data_clean_cases_final.csv")[,-1]
covid_hosp_data <- read.csv("outcome data clean/covid_data_clean_hosp_final.csv")[,-1]
covid_death_data <- read.csv("outcome data clean/covid_data_clean_death_final.csv")[,-1]

#Make dataframes for each outcome (cases, hosp, deaths)
covid_cases <- covid_case_data %>% group_by(age_group) %>% summarise(total_cases = n())
covid_hosp <- covid_hosp_data %>% group_by(age_group) %>% summarise(total_hosp = n()) %>% mutate(total_hosp = total_hosp * (1/0.75))
covid_death <- covid_death_data %>% group_by(age_group) %>% summarise(total_death = n())

#Merge and calculate rates
merged_data <- melt(merge(merge(covid_cases, covid_hosp, by = "age_group", all.x = TRUE),
                     covid_death, by = "age_group", all.x = TRUE) %>%
  mutate(case_fatality = round(total_death/total_cases * 100, 2),
         case_hospitalization = round(total_hosp/total_cases * 100, 2)) %>%
  select(age_group, case_fatality, case_hospitalization))

#Plot
ggplot(merged_data, aes(age_group, value, fill = variable)) +
  geom_col(position = position_dodge()) +
  coord_flip() +
  xlab("Age Group") +
  ylab("Percent (%)") +
  scale_fill_discrete(name = "Rate", labels = c("Case Fatality Rate", "Case Hospitalization Rate"))+
  ggtitle("Case Fatality Rate and Case Hospitalization Rate by Age over Last 6 Months") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))

#########################################################
#Model Calibration/Prediction Figure

#Read more dataframes
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

#Regression Model: Simple Model using only age_group and vaccination status as predictors
calibration_data <- covid_death %>% filter(weeks_since_july2022 %in% c(0:26))

#Model Calibration (on 6 months of data)
model1 <- glm(num_death_adj ~ age_group + boost_2_vax_status_match,
              data = calibration_data, "quasipoisson")
summary(model1)


#Model Calibration Check (Predicting on 6 months of data that it was calibrated on)
validation_period <- covid_death %>% filter(weeks_since_july2022 %in% c(0:26))
pred <- predict(model1, newdata = validation_period, se.fit = T)
validation_period$prediction <- exp(pred$fit)

model_cases <- melt(validation_period %>% group_by(weeks_since_july2022) %>% 
                      summarise(observed = sum(num_cases_adj),prediction = sum(prediction)) %>%
                      mutate(observed = cumsum(observed), prediction = cumsum(prediction)),
                    id = "weeks_since_july2022")

model_hosp <- melt(validation_period %>% group_by(weeks_since_july2022) %>% 
                     summarise(observed = sum(num_hosp_adj),prediction = sum(prediction)) %>%
                     mutate(observed = cumsum(observed), prediction = cumsum(prediction)),
                   id = "weeks_since_july2022")

model_death <- melt(validation_period %>% group_by(weeks_since_july2022) %>% 
                      summarise(observed = sum(num_death_adj),prediction = sum(prediction)) %>%
                      mutate(observed = cumsum(observed), prediction = cumsum(prediction)),
                    id = "weeks_since_july2022")

p1 <- ggplot(model_cases, aes(x = weeks_since_july2022, y = value, color = variable)) +
  geom_line(size = 1.5)+
 # xlab("Number of Weeks Since July 23, 2022") +
 # ylab("Total Outcomes") +
  scale_color_discrete(labels = c("Observed Counts", "Predicted Counts"))+
  scale_y_continuous(label=comma) +
  ggtitle("Cases") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.title= element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank(),
        legend.position = "none")
p1

p2 <- ggplot(model_hosp, aes(x = weeks_since_july2022, y = value, color = variable)) +
  geom_line(size = 1.5)+
  # xlab("Number of Weeks Since July 23, 2022") +
  # ylab("Total Outcomes") +
  scale_color_discrete(labels = c("Observed Counts", "Predicted Counts"))+
  scale_y_continuous(label=comma) +
  ggtitle("Hospitalizations") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.title= element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank(),
        legend.position = "none")
p2

p3 <- ggplot(model_death, aes(x = weeks_since_july2022, y = value, color = variable)) +
  geom_line(size = 1.5)+
  # xlab("Number of Weeks Since July 23, 2022") +
  # ylab("Total Outcomes") +
  scale_color_discrete(labels = c("Observed Counts", "Predicted Counts"))+
  scale_y_continuous(label=comma) +
  ggtitle("Deaths") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.title= element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank())
p3

grid.arrange(
  p1,p2,p3,
  nrow = 1,
  bottom = "Number of Weeks Since July 23, 2022",
  left = "Total Outcomes",
  layout_matrix = rbind(c(1, 1,1,2, 2, 2, 3,3, 3,3)))

#########################################################
#High continued baseline vaccine uptake Figure

uptake_df <- melt(data.frame(weeks = c(0:50),
                        background_vaccination = rep(.44, 51),
                        main_analysis_vaccinated = c(rep(NA, 25), rep(.7, 26))) %>%
  mutate(background_vaccination = case_when((weeks %in% c(25:29)) ~ background_vaccination + 0.033,
                   (weeks %in% c(30:33)) ~ background_vaccination + 0.066,
                   (weeks %in% c(34:38)) ~ background_vaccination + 0.099,
                   (weeks %in% c(39:43)) ~ background_vaccination + 0.1333,
                   (weeks %in% c(44:47)) ~ background_vaccination + 0.166,
                   (weeks %in% c(48:50)) ~ background_vaccination + 0.198,
                   TRUE ~ background_vaccination )), id = "weeks")

ggplot(uptake_df %>% filter(!is.na(value)), aes(x = weeks, y = value * 100, color = variable)) +
  geom_line()+
  xlab("Number of Weeks Since July 23, 2022") +
  ylab("Vaccination Coverage (%)") +
  scale_color_discrete(name = "Analysis Comparison",
                       labels = c("Control with additional \n background vaccination", 
                                  "Primary analysis \n No prospective vaccination"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.spacing.y = unit(0.5, 'cm'))+
  ylim(0, 100) +
  geom_vline(xintercept = 25, linetype = "dashed") +
  annotate(geom = "text",
           label = "Prediction Period Starts",
           x = 23,
           y = 65,
           angle = 90, 
           vjust = 1) +
  guides(color = guide_legend(byrow = TRUE))

#########################################################
#Outcomes over time (Figure 1, Main Manuscript)

#read in data
cases <- read.csv("outcome data clean/covid_data_clean_cases_all.csv")[,-1]
hosp <- read.csv("outcome data clean/covid_data_clean_hosp_all.csv")[,-1]
death <- read.csv("outcome data clean/covid_data_clean_death_all.csv")[,-1]
vaccines <- read.csv('outcome data clean/vaccines_all.csv')[,-1]

#clean data
cases_clean <- cases %>%
  mutate(dtepisode = as.Date(dtepisode),
         weeks_since_jan2020 = as.integer(floor(difftime(strptime(dtepisode, format="%Y-%m-%d"),strptime("2020-01-01", format="%Y-%m-%d"), units = 'weeks')))) %>%
  group_by(weeks_since_jan2020) %>% summarise(total = sum(total_cases))

hosp_clean <- hosp %>%
  mutate(dt_hosp = as.Date(dt_hosp),
         weeks_since_jan2020 = as.integer(floor(difftime(strptime(dt_hosp, format="%Y-%m-%d"),strptime("2020-01-01", format="%Y-%m-%d"), units = 'weeks')))) %>%
  group_by(weeks_since_jan2020) %>% summarise(total = sum(total_hosp))

death_clean <- death %>%
  mutate(dt_death = as.Date(dt_death),
         weeks_since_jan2020 = as.integer(floor(difftime(strptime(dt_death, format="%Y-%m-%d"),strptime("2020-01-01", format="%Y-%m-%d"), units = 'weeks')))) %>%
  group_by(weeks_since_jan2020) %>% summarise(total = sum(total_death))

vaccine_clean <- melt(vaccines %>%
  mutate(administered_date = as.Date(administered_date),
         weeks_since_dec2020 = as.integer(floor(difftime(strptime(administered_date, format="%Y-%m-%d"),strptime("2020-12-01", format="%Y-%m-%d"), units = 'weeks')))) %>%
  group_by(weeks_since_dec2020)  %>% summarise(cumulative_at_least_one_dose = max(cumulative_at_least_one_dose),
                                               cumulative_fully_vaccinated = max(cumulative_fully_vaccinated),
                                               cumulative_booster_recip_count = max(cumulative_booster_recip_count),
                                               cumulative_bivalent_booster_recip_count = max(cumulative_bivalent_booster_recip_count),
                                               booster_eligible = max(booster_eligible),
                                               bivalent_eligible = max(bivalent_eligible)) %>%
  mutate(unvaccinated = round((39240000 - cumulative_at_least_one_dose)/39240000 * 100, 2),
         primary_series = round((cumulative_fully_vaccinated/39240000) * 100, 2),
         booster = round((cumulative_booster_recip_count/booster_eligible) * 100, 2),
         bivalent = round((cumulative_bivalent_booster_recip_count/bivalent_eligible) * 100, 2)), id = "weeks_since_dec2020")

#cases plot
p1 <- ggplot(data = cases_clean, aes(x = weeks_since_jan2020, y = total)) +
  geom_line(size = 1.5, color = 'cyan4') +
  xlab("Time") + ylab("Weekly Cases") +
  ggtitle("A. Cases")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text=element_text(family="Times New Roman", face="bold", size = 12),
        axis.text.x = element_text(angle=90, hjust=1)) +
  scale_y_continuous(labels=comma) +
  scale_x_continuous(breaks=c(0, 13, 30, 48, 61, 74, 91, 109, 126, 143, 156),
                     labels = c("Jan 2020","Apr 2020","Aug 2020","Dec 2020","Mar 2021", "Jun 2021",
                                "Nov 2021","Feb 2021","Jun 2022", "Oct 2022","Jan 2023")) +
  annotate("rect", xmin = 133, xmax = 159, ymin = 0, ymax = Inf,
           alpha = .1, fill = 'blue')

p1

#hosp plot
p2 <- ggplot(data = hosp_clean, aes(x = weeks_since_jan2020, y = total)) +
  geom_line(size = 1.5, color = 'cyan4') +
  xlab("Time") + ylab("Weekly Hospitalizations") +
  ggtitle("B. Hospitalizations")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text=element_text(family="Times New Roman", face="bold", size = 12),
        axis.text.x = element_text(angle=90, hjust=1)) +
  scale_y_continuous(labels=comma) +
  scale_x_continuous(breaks=c(0, 13, 30, 48, 61, 74, 91, 109, 126, 143, 156),
                     labels = c("Jan 2020","Apr 2020","Aug 2020","Dec 2020","Mar 2021", "Jun 2021",
                                "Nov 2021","Feb 2021","Jun 2022", "Oct 2022","Jan 2023"))+
  annotate("rect", xmin = 133, xmax = 159, ymin = 0, ymax = Inf,
           alpha = .1, fill = 'blue')


p2

#death plot
p3 <- ggplot(data = death_clean, aes(x = weeks_since_jan2020, y = total)) +
  geom_line(size = 1.5, color = 'cyan4') +
  xlab("Time") + ylab("Weekly Deaths") +
  ggtitle("C. Deaths")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text=element_text(family="Times New Roman", face="bold", size = 12),
        axis.text.x = element_text(angle=90, hjust=1)) +
  scale_y_continuous(labels=comma) +
  scale_x_continuous(breaks=c(0, 13, 30, 48, 61, 74, 91, 109, 126, 143, 156),
                     labels = c("Jan 2020","Apr 2020","Aug 2020","Dec 2020","Mar 2021", "Jun 2021",
                                "Nov 2021","Feb 2021","Jun 2022", "Oct 2022","Jan 2023"))+
  annotate("rect", xmin = 133, xmax = 159, ymin = 0, ymax = Inf,
           alpha = .1, fill = 'blue')


p3

#vaccination plot
p4 <- ggplot(data = vaccine_clean %>% filter(variable %in% c("unvaccinated", "primary_series", "booster", "bivalent")), aes(x = weeks_since_dec2020, y = value, color = variable)) +
  geom_line(size = 1.5)  +
  xlab("Time") + ylab("Vaccination Coverage (%)") +
  ggtitle("D. Vaccination")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text=element_text(family="Times New Roman", face="bold", size = 12),
        axis.text.x = element_text(angle=90, hjust=1)) +
  scale_y_continuous(labels=comma) +
  scale_x_continuous(breaks=c(0, 9, 17, 26, 35, 43, 52, 61, 69, 78, 87, 95, 104, 111),
                     labels = c("Dec 2020","Feb 2021","Apr 2021", "Jun 2021", "Aug 2021",
                                "Oct 2021", "Dec 2021", "Feb 2022", "Apr 2022", "Jun 2022",
                                "Aug 2022", "Oct 2022", "Dec 2022", "Feb 2023"))+
  annotate("rect", xmin = 85, xmax = 111, ymin = 0, ymax = Inf,
           alpha = .1, fill = 'blue') +
  scale_color_discrete(name = "Vaccination Status",
                       labels = c("Unvaccinated", "Fully Vaccinated", "Boosted", "Boosted (Bivalent)"))


p4

grid.arrange(
  p1,p2,p3,p4,
  nrow = 2, layout_matrix = rbind(c(1, 2, 3),
                                  c(4, 4, 4)))

#########################################################
#Comparing reduction in Effectiveness in NNT (Figure S1, Appendix)

cases_vax <- read.csv("sensitivity-reduction in effectiveness/data/results-vax-cases.csv")[,-1]
hosp_vax <- read.csv("sensitivity-reduction in effectiveness/data/results-vax-hosp.csv")[,-1]
death_vax <- read.csv("sensitivity-reduction in effectiveness/data/results-vax-death.csv")[,-1]

hosp_pax <- read.csv("sensitivity-reduction in effectiveness/data/results-pax-hosp.csv")[,-1]
death_pax <- read.csv("sensitivity-reduction in effectiveness/data/results-pax-death.csv")[,-1]

combined_vax <- melt(data.frame(perc_reduction = c(1:1000)/10,
                       cases = cases_vax$NNT_diff,
                       hosp = hosp_vax$NNT_diff,
                       death = death_vax$NNT_diff), id = "perc_reduction")

combined_pax <- melt(data.frame(perc_reduction = c(1:1000)/10,
                                hosp = hosp_pax$NNT_diff,
                                death = death_pax$NNT_diff), id = "perc_reduction")


p1 <- ggplot(combined_vax, aes(x = perc_reduction, y = value, color = variable)) +
  geom_line(size = 1.5)+
  xlab("Percent Reduction in Vaccine Effectiveness (%)") +
  ylab("NNT Difference") +
  scale_color_discrete(labels = c("Cases", "Hospitalizations", "Deaths"), name = "Outcome")+
  scale_y_continuous(label=comma, limits = c(-50000, 17000)) +
  ggtitle("A. Bivalent Vaccines") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.title= element_blank()) +
  geom_hline(yintercept = 0, linetype = "dashed")

p1

p2 <- ggplot(combined_pax, aes(x = perc_reduction, y = value, color = variable)) +
  geom_line(size = 1.5)+
  xlab("Percent Reduction in Paxlovid Effectiveness (%)") +
  ylab("NNT Difference") +
  scale_color_discrete(labels = c("Hospitalizations", "Deaths"), name = "Outcome")+
  scale_y_continuous(label=comma, limits = c(-400, 40)) +
  ggtitle("B. Nirmatrelvir-Ritonavir Treatment") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.title= element_blank()) +
  geom_hline(yintercept = 0, linetype = "dashed")

p2


grid.arrange(
  p1,p2,
  nrow = 1)
