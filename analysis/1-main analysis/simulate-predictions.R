###################################################################################################
#Title: Parameter Simulation
#Author: Hailey Park
#Date: January 10, 2023
#Taken from Sophia's code
###################################################################################################

rm(list=ls())

setwd("/mnt/projects/covid_partners/ucsf_lo/Bivalent Booster Paxlovid Impact Analysis/final-code")


#Define functions
objective.function <- function(params) {
  alpha <- params[1]
  beta <- params[2]
  calculated.quantiles <- qbeta(p=c(0.025, 0.975), shape1=alpha, shape2=beta)
  squared.error.quantiles <- sum((intended.quantiles - calculated.quantiles)^2)
  
  calculated.mean <- calculate.mean(alpha, beta)
  squared.error.mean <- (intended.mean - calculated.mean)^2
  
  return(squared.error.quantiles + squared.error.mean)
}


calculate.mean <- function(alpha, beta) {
  return(alpha/(alpha+beta))
}


optimal <- function() {
  starting.params <- c(5, 1)
  nlm.result <- nlm(f = objective.function, p = starting.params)
  optimal.alpha <- nlm.result$estimate[1]
  optimal.beta <- nlm.result$estimate[2]
  print(c(alpha=optimal.alpha, beta=optimal.beta))
  c(optimal.alpha, optimal.beta)
}

# variable parameters
# fix n as number of simulations
n <- 1000

ve_parameters <- data.frame(trial = 1:n)
pax_parameters <- data.frame(trial = 1:n)
ve_waning_parameters <- data.frame(trial = 1:n)

#Run for each VE and Paxlovid effectiveness estimate
set.seed(42)

p <- matrix(runif(n*10, 0, 1), nrow = n, ncol=10)

  #VE for Unvaccinated against Cases
intended.quantiles <- c(0.29, 0.35)
intended.mean <- 0.32
opt_unvax_case <- optimal()
ve_parameters$unvax_case <- qbeta(p[,1], opt_unvax_case[1], opt_unvax_case[2])

  #VE for Unvaccinated against Hosp
intended.quantiles <- c(0.21, 0.34)
intended.mean <- 0.29
opt_unvax_hosp <- optimal()
ve_parameters$unvax_hosp <- qbeta(p[,2], opt_unvax_hosp[1], opt_unvax_hosp[2])

  #VE for Unvaccinated against Death
intended.quantiles <- c(0.46, 0.69)
intended.mean <- 0.57
opt_unvax_death <- optimal()
ve_parameters$unvax_death <- qbeta(p[,3], opt_unvax_death[1], opt_unvax_death[2])

  #VE for Primary series against Cases
intended.quantiles <- c(0.25, 0.53)
intended.mean <- 0.41
opt_prim_case <- optimal()
ve_parameters$prim_case <- qbeta(p[,4], opt_prim_case[1], opt_prim_case[2])

  #VE for Primary series against Hosp
intended.quantiles <- c(0.41, 0.69)
intended.mean <- 0.57
opt_prim_hosp <- optimal()
ve_parameters$prim_hosp <- qbeta(p[,5], opt_prim_hosp[1], opt_prim_hosp[2])

  #VE for Primary series against Death
intended.quantiles <- c(0.67, 0.83)
intended.mean <- 0.75
opt_prim_death <- optimal()
ve_parameters$prim_death <- qbeta(p[,6], opt_prim_death[1], opt_prim_death[2])

  #VE for Boosted 1 dose against Cases
intended.quantiles <- c(0.01, 0.44)
intended.mean <- 0.26
opt_boost_case <- optimal()
ve_parameters$boost_case <- qbeta(p[,7], opt_boost_case[1], opt_boost_case[2])

  #VE for Boosted 1 dose against Hosp
intended.quantiles <- c(0.13, 0.56)
intended.mean <- 0.38
opt_boost_hosp <- optimal()
ve_parameters$boost_hosp <- qbeta(p[,8], opt_boost_hosp[1], opt_boost_hosp[2])

  #VE for Boosted 1 dose against Death
intended.quantiles <- c(0.34, 0.90)
intended.mean <- 0.62
opt_boost_death <- optimal()
ve_parameters$boost_death <- qbeta(p[,9], opt_boost_death[1], opt_boost_death[2])

#VE for Boosted 2 dose against Cases
intended.quantiles <- c(0.32, 0.47)
intended.mean <- 0.40
opt_boost_case <- optimal()
ve_parameters$boost2_case <- qbeta(p[,7], opt_boost_case[1], opt_boost_case[2])

#VE for Boosted 2 dose against Hosp
intended.quantiles <- c(0.12, 0.78)
intended.mean <- 0.56
opt_boost_hosp <- optimal()
ve_parameters$boost2_hosp <- qbeta(p[,8], opt_boost_hosp[1], opt_boost_hosp[2])

#VE for Boosted 2 dose against Death
intended.quantiles <- c(0.27, 0.81)
intended.mean <- 0.63
opt_boost_death <- optimal()
ve_parameters$boost2_death <- qbeta(p[,9], opt_boost_death[1], opt_boost_death[2])


  #Pax for Unvaccinated against Hosp
intended.quantiles <- c(0.51, 1)
intended.mean <- 0.89
opt_unvax_hosp <- optimal()
pax_parameters$unvax_hosp <- qbeta(p[,1],opt_unvax_hosp[1], opt_unvax_hosp[2])

  #Pax for Unvaccinated against Death
intended.quantiles <- c(0.51, 1)
intended.mean <- 0.89
opt_unvax_death <- optimal()
pax_parameters$unvax_death <- qbeta(p[,2], opt_unvax_death[1], opt_unvax_death[2])


  #Pax for Vaccinated against Hosp
intended.quantiles <- c(0.19, 0.56)
intended.mean <- 0.40
opt_vax_hosp <- optimal()
pax_parameters$vax_hosp <- qbeta(p[,3], opt_vax_hosp[1], opt_vax_hosp[2])

  #Pax for Vaccinated against Death
intended.quantiles <- c(0.29, 0.88)
intended.mean <- 0.71
opt_vax_death <- optimal()
pax_parameters$vax_death <- qbeta(p[,4], opt_vax_death[1], opt_vax_death[2])

  #Waning for Unvax 0-4 weeks
intended.quantiles <- c(0.417, 0.487)
intended.mean <- 0.4535
opt_unvax_0_to_4 <- optimal()
ve_waning_parameters$unvax_0_to_4 <- qbeta(p[,1], opt_unvax_0_to_4[1], opt_unvax_0_to_4[2])

  #Waning for Unvax >5 weeks
intended.quantiles <- c(0.286, 0.346)
intended.mean <- 0.317
opt_unvax_5_plus <- optimal()
ve_waning_parameters$unvax_5_plus <- qbeta(p[,2], opt_unvax_5_plus[1], opt_unvax_5_plus[2])

  #Waning for Primary Series 0-4 weeks
intended.quantiles <- c(0.514, 0.598)
intended.mean <- 0.558
opt_prim_0_to_4 <- optimal()
ve_waning_parameters$prim_0_to_4 <- qbeta(p[,3], opt_prim_0_to_4[1], opt_prim_0_to_4[2])

  #Waning for Primary Series 5-9 weeks
intended.quantiles <- c(0.389, 0.465)
intended.mean <- 0.428
opt_prim_5_to_9 <- optimal()
ve_waning_parameters$prim_5_to_9 <- qbeta(p[,4], opt_prim_5_to_9[1], opt_prim_5_to_9[2])

  #Waning for Primary Series 10-13 weeks
intended.quantiles <- c(0.227, 0.286)
intended.mean <- 0.257
opt_prim_10_to_13 <- optimal()
ve_waning_parameters$prim_10_to_13 <- qbeta(p[,5], opt_prim_10_to_13[1], opt_prim_10_to_13[2])

  #Waning for Primary Series 14-17 weeks
intended.quantiles <- c(0.173, 0.232)
intended.mean <- 0.203
opt_prim_14_to_17 <- optimal()
ve_waning_parameters$prim_14_to_17 <- qbeta(p[,6], opt_prim_14_to_17[1], opt_prim_14_to_17[2])

  #Waning for Primary Series >18 weeks
intended.quantiles <- c(0.158, 0.224)
intended.mean <- 0.192
opt_prim_18_plus <- optimal()
ve_waning_parameters$prim_18_plus <- qbeta(p[,7], opt_prim_18_plus[1], opt_prim_18_plus[2])

  #Waning for Boosted 0-5 weeks
intended.quantiles <- c(0.23, 0.51)
intended.mean <- 0.39
opt_boost_0_to_5 <- optimal()
ve_waning_parameters$boost_0_to_5 <- qbeta(p[,8], opt_boost_0_to_5[1], opt_boost_0_to_5[2])

  #Waning for Boosted 6-15 weeks
intended.quantiles <- c(0.16, 0.45)
intended.mean <- 0.32
opt_boost_6_to_15 <- optimal()
ve_waning_parameters$boost_6_to_15 <- qbeta(p[,9], opt_boost_6_to_15[1], opt_boost_6_to_15[2])

  #Waning for Boosted 16plus weeks
intended.quantiles <- c(0, 0.17)
intended.mean <- 0
opt_boost_16_plus <- optimal()
ve_waning_parameters$boost_16_plus <- qbeta(p[,10], opt_boost_16_plus[1], opt_boost_16_plus[2])



inspection <- ve_parameters
inspection %>% apply(2, quantile, probs=c(0.025, 0.975))

saveRDS(ve_parameters, "data/simulated-ve-parameters.RDS")
saveRDS(pax_parameters, "data/simulated-pax-parameters-extra.RDS")
saveRDS(ve_waning_parameters, "data/simulated-ve-waning-parameters.RDS")

