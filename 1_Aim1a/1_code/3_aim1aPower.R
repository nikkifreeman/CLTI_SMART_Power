# 1_aim1aPower.R
# Author: Nikki Freeman
# Last modified: 25 April 2022

# Preamble ---------------------------------------------------------------------
# Code to generate the simulation analysis for Aim 1a
#
# Required packages: tidyverse, tmle, SuperLearner, tictoc

# Packages and scripts ---------------------------------------------------------
library(tidyverse)
library(tictoc)
library(tmle)
library(SuperLearner)

if(str_detect(getwd(), "resubmission")){
  source("./1_Aim1a/1_code/1_generateSampleAim1a.R")
  source("./1_Aim1a/1_code/0_universalParameters.R")
} else{
  source("1_generateSampleAim1a.R")
  source("0_universalParameters.R")
}

# Calculate the power for different values of the parameters -------------------
# vectors to hold on to output
n_out <- p_out <- l_out <- psi_out <- var_psi_out <- CI_lower_out <- CI_upper_out <- pvalue_out <- vector(length = 0)

set.seed(10222019) # Muffie's birthday
# n keeps track of the sample size
for(n in N){
  # p keeps track of the parameter values
  for(p in 1:nrow(parameterGrid)){
    print(c(n,p))
    # l keeps track of the number of simulated data sets
    for(l in 1:L){
      # Generate a data set
      simData <- simulateDataAim1a_exponential(
        N = n, p_trt1 = 0.5, 
        mort0 = parameterGrid$mort2yr0[p], 
        mort1 = parameterGrid$mort2yr1[p], 
        recur0 = parameterGrid$woundRecur0[p], 
        recur1 = parameterGrid$woundRecur1[p], 
        C_L = C_L, 
        beta_mort = rep(0, 4), 
        beta_recur = rep(0, 4)) %>%
        mutate(observed = 1 - ltfu)
      
      # Calculate ATE
      tmle_out <- tmle(
        Y = simData$woundFreeAliveDays,
        A = simData$A1,
        W = simData[ ,c("X1", "X2", "X3", "X4")],
        Delta = simData$observed)
      
      temp <- summary(tmle_out)$estimates$ATE
      n_out <- c(n_out, n)
      p_out <- c(p_out, p) 
      l_out <- c(l_out, l)
      psi_out <- c(psi_out, temp$psi)
      var_psi_out <- c(var_psi_out, temp$var.psi)
      CI_lower_out <- c(CI_lower_out, temp$CI[1]) 
      CI_upper_out <- c(CI_upper_out, temp$CI[2])
      pvalue_out <- c(pvalue_out, temp$pvalue)
    }
  }
}

out <- data.frame(n = n_out, p = p_out, l = l_out,
                  psi = psi_out, var_psi = var_psi_out,
                  CI_lower = CI_lower_out, CI_upper = CI_upper_out,
                  pvalue = pvalue_out)
readr::write_csv(x = out, file = "../2_pipeline/3_aim1aPower.csv")