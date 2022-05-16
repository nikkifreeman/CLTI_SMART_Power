# 0_generateSample.R
# Author: Nikki Freeman
# Last modified: 25 April 2022

# Preamble ---------------------------------------------------------------------
# These are the functions for generating a single sample for the simulation
# analysis. 
# Required packages: tidyverse

#' Generate baseline data
#'
#' @param N Positive integer; corresponds to the size of the sample
#'
#' @return dataframe; columns are id, X1 ~ N(0, 1), X2 ~ N(0, 1), 
#' X3 ~ Bernoulli(0.6), and X4 ~ Bernoulli(0.4)
generateBaselineData <- function(N){
  # Patient ids
  id <- 1:N
  
  # Covariates
  X1 <- rnorm(n = N)
  X2 <- rnorm(n = N)
  X3 <- rbinom(n = N, size = 1, prob = 0.6)
  X4 <- rbinom(n = N, size = 1, prob = 0.4)
  data.frame(id, X1, X2, X3, X4)
}

#' Assign first randomization 
#'
#' @param df dataframe; contains the baseline data
#' @param N positive integer; should be the same N used to generate baseline data
#' @param p_trt1 real in (0, 1); probability of recieving "treatment 1"
#'
#' @return dataframe; baseline data with the randomized first treatment
assignFirstRandomization <- function(df, N, p_trt1){
  trt1 <- sample(1:N, size = N*p_trt1, replace = FALSE)
  df %>% 
    mutate(A1 = if_else(id %in% trt1, 1, 0))
}

#' Calculate Lambda (exponential distribution)
#'
#' Helper function to calculate the lambda for the exponential distribution
#' @param S real in (0,1); probability of survival over t
#' @param t integer; time over which S is measured
#'
#' @return real
calculateLambda_exponential <- function(S, t){
  H_0 <- -1*log(1-S)
  H_0/t
}

#' Generate exponentially distributed event time
#'
#' @param S real in (0,1); probability of survival over t
#' @param t integer; time over which S is measured
#' @param expLinPred rea; linear predictor (covariates * beta)
#'
#' @return real; event time
generateEventTime_exponential <- function(S, t, expLinPred){
  lambda <- calculateLambda_exponential(S, t)
  -1*log(runif(1))/(lambda * expLinPred)
}

#' Generate a single simulation data set for Aim 1.1
#'
#' @param N positive integer; same N used for baseline data, sample size
#' @param p_trt1 real in (0, 1); probability of recieving "treatment 1"
#' @param mort0 real in (0, 1); Pr(2 year mortality) in arm 0
#' @param mort1 real in (0, 1); Pr(2 year mortality) in arm 1
#' @param recur0 real in (0, 1); Pr(recurrence in 2 years) in arm 0
#' @param recur1 real in (0, 1); Pr(recurrence in 2 years) in arm 1
#' @param C_L real in (0, 1); Pr(ltfu over the study period)
#' @param beta_mort real[4]; coefficients for X1:X4 in the TTE mortality 
#' @param beta_recur real[4]; coefficients for X1:X4 in the TTE recurrence
#'
#' @return dataframe; contains 1 simulation data set for aim 1.1
simulateDataAim1a_exponential <- function(N, p_trt1, mort0, mort1, recur0, recur1, C_L, beta_mort, beta_recur){
  ltfu_ids <- sample(1:N, size = round(N*C_L), replace = FALSE)
  generateBaselineData(N) %>%
    assignFirstRandomization(., N = N, p_trt1 = p_trt1) %>%
    mutate(linPred_mort = beta_mort[1]*X1 + beta_mort[2]*X2 + beta_mort[3]*X3 + beta_mort[4]*X4,
           expLinPred_mort = exp(linPred_mort)) %>%
    mutate(linPred_recur = beta_recur[1]*X1 + beta_recur[2]*X2 + beta_recur[3]*X3 + beta_recur[4]*X4,
           expLinPred_recur = exp(linPred_recur)) %>%
    rowwise() %>%
    mutate(timeToDeath = 
             if_else(A1 == 0, 
                     generateEventTime_exponential(mort0, t = 365*2, expLinPred_mort), 
                     generateEventTime_exponential(mort1, t = 365*2, expLinPred_mort))) %>%
    mutate(timeToRecur = 
             if_else(A1 == 0,
                     generateEventTime_exponential(recur0, t = 365*2, expLinPred_recur),
                     generateEventTime_exponential(recur1, t = 365*2, expLinPred_recur))) %>%
    mutate(timeToHealFirstWound = truncnorm::rtruncnorm(1, a = 0.5*30, mean = 4.5*30, sd = 2*30)) %>%
    mutate(timeToHealSecondWound = truncnorm::rtruncnorm(1, a = 0.5*30, mean = 4.5*30, sd = 2*30)) %>%
    ungroup() %>%
    mutate(ltfu = if_else(id %in% ltfu_ids, 1, 0)) %>%
    mutate(firstWoundHealingDay = timeToHealFirstWound) %>% 
    rowwise() %>%
    mutate(timeToRecur = max(timeToRecur, firstWoundHealingDay + 1)) %>%
    ungroup() %>%
    mutate(secondWoundHealingDay = timeToRecur + timeToHealSecondWound) %>%
    mutate(outcome_qual = if_else(timeToDeath <= 365*2 & 
                                    timeToDeath <= firstWoundHealingDay, 
                                  "died before healing first wound", "?"),
           outcome_qual = if_else(outcome_qual == "?" & timeToDeath <= 365*2 & 
                                    timeToDeath <= timeToRecur, 
                                  "healed first wound, died before recurrence",
                                  outcome_qual),
           outcome_qual = if_else(outcome_qual == "?" & timeToDeath > 365*2 & 
                                    timeToRecur > 365*2 & 
                                    firstWoundHealingDay <= 365*2, 
                                  "healed first wound, no recurrence, alive", 
                                  outcome_qual),
           outcome_qual = if_else(outcome_qual == "?" & timeToDeath <= 365*2 & 
                                    timeToDeath > firstWoundHealingDay & 
                                    timeToDeath > timeToRecur & 
                                    timeToDeath <= secondWoundHealingDay, 
                                  "healed first wound, recurred, died before healing second wound", 
                                  outcome_qual),
           outcome_qual = if_else(outcome_qual == "?" & timeToDeath <= 365*2 & 
                                    timeToDeath > firstWoundHealingDay & 
                                    timeToDeath > timeToRecur & 
                                    timeToDeath > secondWoundHealingDay, 
                                  "healed first wound, recurred, healed second wound, died", 
                                  outcome_qual),
           outcome_qual = if_else(outcome_qual == "?" & timeToDeath > 365*2 & 
                                    timeToDeath > firstWoundHealingDay & 
                                    timeToDeath > timeToRecur & 
                                    secondWoundHealingDay > 365*2, 
                                  "healed first, recurred, didn't heal second",  
                                  outcome_qual),
           outcome_qual = if_else(outcome_qual == "?" & timeToDeath > 365*2 & 
                                    timeToDeath > firstWoundHealingDay & 
                                    timeToDeath > timeToRecur & 
                                    secondWoundHealingDay <= 365*2, 
                                  "healed first, recurred, healed second", 
                                  outcome_qual))  %>%
    mutate(woundFreeAliveDays = if_else(outcome_qual == "died before healing first wound", 0, -100),
           woundFreeAliveDays = if_else(outcome_qual == "healed first wound, died before recurrence", 
                                        timeToDeath - firstWoundHealingDay, woundFreeAliveDays),
           woundFreeAliveDays = if_else(outcome_qual == "healed first wound, no recurrence, alive", 
                                        365*2 - firstWoundHealingDay, woundFreeAliveDays),
           woundFreeAliveDays = if_else(outcome_qual == "healed first wound, recurred, died before healing second wound", 
                                        timeToRecur - firstWoundHealingDay, woundFreeAliveDays),
           woundFreeAliveDays = if_else(outcome_qual == "healed first wound, recurred, healed second wound, died", 
                                        timeToDeath - timeToHealFirstWound - timeToHealSecondWound, woundFreeAliveDays),
           woundFreeAliveDays = if_else(outcome_qual == "healed first, recurred, didn't heal second", 
                                        timeToRecur - firstWoundHealingDay, woundFreeAliveDays),
           woundFreeAliveDays = if_else(outcome_qual == "healed first, recurred, healed second", 
                                        365*2 - timeToHealFirstWound - timeToHealSecondWound, woundFreeAliveDays))
  
}