# Aim 1a with amputations power
# Author: Nikki Freeman
# Last modified: 23 June 2022

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
  source("./1_Aim1a/1_code/4_generateSampleWithAmpAim1a.R")
  simulationParameters <- read_csv("./1_Aim1a/1_code/aim1a_simulationParameters.csv", show_col_types = FALSE)
} else{
  source("1_generateSampleAim1a.R")
  source("0_universalParameters.R")
  source("4_generateSampleWithAmpAim1a.R")
  simulationParameters <- read_csv("aim1a_simulationParameters.csv")
}

# Command line arguments -------------------------------------------------------
command_args <- commandArgs(trailingOnly = TRUE)
command_args <- as.numeric(command_args)
studySize <- command_args[1]
p <- command_args[2]

# Calculate the power for different values of the parameters -------------------
# df to hold on to output
out <- data.frame()
tic()
set.seed(10222019) # Muffie's birthday
# n keeps track of the sample size
# for(n in command_args[1]){
  # # p keeps track of the parameter values
  # for(p in 1:nrow(simulationParameters)){
    # l keeps track of the number of simulated data sets
    for(l in 1:L){
      # Generate a data set
      simData <- simulateDataAim1aWithAmp_exponential(
        N = studySize, p_trt1 = 0.5, 
        mort0 = simulationParameters$mort2yr0[p], 
        mort1 = simulationParameters$mort2yr1[p], 
        recur0 = simulationParameters$woundRecur0[p], 
        recur1 = simulationParameters$woundRecur1[p], 
        amp0 = simulationParameters$amp0[p],
        amp1 = simulationParameters$amp1[p],
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
      out <- bind_rows(out, 
                       data.frame(N = studySize,
                                  l = l, psi = temp$psi, var_psi = temp$var.psi,
                                  CI_lower = temp$CI[1], CI_upper = temp$CI[2],
                                  pvalue = temp$pvalue,   
                                  mort0 = simulationParameters$mort2yr0[p], 
                                  mort1 = simulationParameters$mort2yr1[p], 
                                  recur0 = simulationParameters$woundRecur0[p], 
                                  recur1 = simulationParameters$woundRecur1[p], 
                                  amp0 = simulationParameters$amp0[p],
                                  amp1 = simulationParameters$amp1[p]))
    }
  # }
# }
toc()



if(str_detect(getwd(), "resubmission")){
  outFileName <- paste0("./1_Aim1a/2_pipeline/8_aim1aWithAmpPower_", studySize, "_", p, ".csv")
} else{
  outFileName <- paste0("../2_pipeline/8_aim1aWithAmpPower_", studySize, "_", p, ".csv")
}
readr::write_csv(x = out, file = outFileName)
