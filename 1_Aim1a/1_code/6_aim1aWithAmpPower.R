# Aim 1a with amputations power
# Author: Nikki Freeman
# Last modified: 9 June 2022

# Preamble ---------------------------------------------------------------------
# Code to generate the simulation analysis for Aim 1a
#
# Required packages: tidyverse, tmle, SuperLearner, tictoc

# Command line arguments -------------------------------------------------------
command_args <- commandArgs(trailingOnly = TRUE)
command_args <- as.numeric(command_args)

# Packages and scripts ---------------------------------------------------------
library(tidyverse)
library(tictoc)
library(tmle)
library(SuperLearner)

if(str_detect(getwd(), "resubmission")){
  source("./1_Aim1a/1_code/1_generateSampleAim1a.R")
  source("./1_Aim1a/1_code/0_universalParameters.R")
  source("./1_Aim1a/1_code/4_generateSampleWithAmpAim1a.R")
} else{
  source("1_generateSampleAim1a.R")
  source("0_universalParameters.R")
  source("4_generateSampleWithAmpAim1a.R")
}

# Calculate the power for different values of the parameters -------------------
# vectors to hold on to output
n_out <- p_out <- l_out <- psi_out <- var_psi_out <- CI_lower_out <- CI_upper_out <- pvalue_out <- vector(length = 0)

tic()
set.seed(10222019) # Muffie's birthday
# n keeps track of the sample size
# for(n in command_args[1]){
  # # p keeps track of the parameter values
  # for(p in 1:nrow(parameterGrid)){
    # l keeps track of the number of simulated data sets
    for(l in 1:L){
      # Generate a data set
      simData <- simulateDataAim1aWithAmp_exponential(
        N = command_args[1], p_trt1 = 0.5, 
        mort0 = parameterGrid$mort2yr0[command_args[4]], 
        mort1 = parameterGrid$mort2yr1[command_args[4]], 
        recur0 = parameterGrid$woundRecur0[command_args[4]], 
        recur1 = parameterGrid$woundRecur1[command_args[4]], 
        amp0 = command_args[2],
        amp1 = command_args[3],
        C_L = C_L, 
        beta_mort = rep(0, 4), 
        beta_recur = rep(0, 4),
        beta_amp = rep(0, 4)) %>%
        mutate(observed = 1 - ltfu)
      
      # Calculate ATE
      tmle_out <- tmle(
        Y = simData$woundFreeAliveDays,
        A = simData$A1,
        W = simData[ ,c("X1", "X2", "X3", "X4")],
        Delta = simData$observed)
      
      temp <- summary(tmle_out)$estimates$ATE
      n_out <- c(n_out, n)
      l_out <- c(l_out, l)
      psi_out <- c(psi_out, temp$psi)
      var_psi_out <- c(var_psi_out, temp$var.psi)
      CI_lower_out <- c(CI_lower_out, temp$CI[1]) 
      CI_upper_out <- c(CI_upper_out, temp$CI[2])
      pvalue_out <- c(pvalue_out, temp$pvalue)
    }
  # }
# }
toc()

out <- data.frame(l = l_out, psi = psi_out, var_psi = var_psi_out,
                  CI_lower = CI_lower_out, CI_upper = CI_upper_out,
                  pvalue = pvalue_out)
out$n <- command_args[1]
out$mort2yr0 <- parameterGrid$mort2yr0[command_args[4]]
out$mort2yr1 <- parameterGrid$mort2yr1[command_args[4]]
out$woundRecur0 <- parameterGrid$woundRecur0[command_args[4]]
out$woundRecur1 <- parameterGrid$woundRecur1[command_args[4]]
out$C_L <- parameterGrid$C_L[command_args[4]]

if(str_detect(getwd(), "resubmission")){
  outFileName <- paste0("./1_Aim1a/2_pipeline/6_aim1aWithAmpPower_n", command_args[1], 
                        "_", round(command_args[2]*100), "_", round(command_args[3]*100),
                                                                    "_", command_args[4], ".csv")
} else{
  outFileName <- paste0("../2_pipeline/6_aim1aWithAmpPower_n", command_args[1], 
                        "_", round(command_args[2]*100), "_", round(command_args[3]*100), "_", command_args[4], ".csv")
}
readr::write_csv(x = out, file = outFileName)
