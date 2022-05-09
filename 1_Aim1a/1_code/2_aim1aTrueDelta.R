# 1_aim1aTrue.R
# Author: Nikki Freeman
# Last modified: 25 April 2022

# Preamble ---------------------------------------------------------------------
# Code to generate a very large number of simulated participants under varying 
# assumptions about the true 2-year mortality and true 2-year wound recurrence
# rates. This can be used to calculate the true difference in the number of
# wound free days between the two different arms of the embedded RCT.
#
# Required packages: tidyverse, tictoc


# Load packages and scripts ----------------------------------------------------
# Packages
library(tidyverse)
library(tictoc)
# scripts
if(str_detect(getwd(), "resubmission")){
  source("./1_Aim1a/1_code/1_generateSampleAim1a.R")
  source("./1_Aim1a/1_code/0_universalParameters.R")
} else{
  source("1_generateSampleAim1a.R")
  source("0_universalParameters.R")
}

# Calculate the true values under various parameters ---------------------------
parameterGrid <- parameterGrid %>%
  add_column(Delta = -10000)

# Generate a large sample (N = 50000) using the parameters and calculate 
# the true difference in the means as an approximation of the true population 
# difference between the two arms (by LLN)

tic()
set.seed(5678)
for(i in 1:nrow(parameterGrid)){
  simData <- simulateDataAim1a_exponential(
    N = 50000, p_trt1 = 0.5, 
    mort0 = parameterGrid$mort2yr0[i], 
    mort1 = parameterGrid$mort2yr1[i], 
    recur0 = parameterGrid$woundRecur0[i], 
    recur1 = parameterGrid$woundRecur1[i], 
    C_L = C_L, 
    beta_mort = rep(0, 4), 
    beta_recur = rep(0, 4)) 
  temp <- simData %>% group_by(A1) %>% summarise(groupMean = mean(woundFreeAliveDays)) 
  parameterGrid$Delta[i] <- temp$groupMean[temp$A1 == 1] - temp$groupMean[temp$A1 == 0]
}
toc()

readr::write_csv(x = parameterGrid, file = "../2_pipeline/2_aim1aTrueDelta.csv")


