# Power analysis for aim 1b
# Author: Nikki Freeman
# Last modified: 5 May 2022
# Dependencies (packages): tidyverse, sandwich
# Dependencies (scripts): 2_generateSampleAim1b.R

# Read the commnand line arguments
command_args <- commandArgs(trailingOnly = TRUE)
command_args <- as.numeric(command_args)
command_args <- c(650, 1, .2, .2, 1)
# Packages ---------------------------------------------------------------------
library(tidyverse)
library(sandwich)
library(tictoc)

# Scripts ----------------------------------------------------------------------
if(str_detect(getwd(), "resubmission")){
  source("./2_Aim1b/1_code/1_aim1bPower_fnx.R")
  source("./2_Aim1b/1_code/2_generateSampleAim1b.R")
  source("./1_Aim1a/1_code/0_universalParameters.R")
} else{
  source("1_aim1bPower_fnx.R")
  source("2_generateSampleAim1b.R")
  source("../../1_Aim1a/1_code/0_universalParameters.R")
}

# Simulation parameters --------------------------------------------------------
studySize <- as.numeric(command_args[1])
designProbs <- rep(0.5, 5)
dominantRegime <- as.numeric(command_args[2])

# Dataframe for output 
out <- data.frame()

# Simulation -------------------------------------------------------------------
tic()
# for(i in 1:nrow(parameterGrid)){
mort2yr_dom <- parameterGrid[command_args[5], "mort2yr1"]
recur2yr_dom <- parameterGrid[command_args[5], "woundRecur1"]
mort2yr_notDom <- parameterGrid[command_args[5], "mort2yr0"]
recur2yr_notDom <- parameterGrid[command_args[5], "woundRecur0"]
amp2yr_dom <- as.numeric(command_args[4])
amp2yr_notDom <- as.numeric(command_args[3])
censorRate <- C_L
currParamDf <- data.frame(mort2yr_dom = mort2yr_dom, 
                          recur2yr_dom = recur2yr_dom,
                          mort2yr_notDom = mort2yr_notDom, 
                          recur2yr_notDom = recur2yr_notDom,
                          amp2yr_dom = amp2yr_dom,
                          amp2yr_notDom = amp2yr_notDom,
                          censorRate = C_L)
for(l in 1:L){
  resultsOneRun <- doOneSimRun(studySize = studySize, 
                               designProbs = designProbs, 
                               responseProbRange = responseProbRange,
                               dominantRegime = dominantRegime, 
                               mort2yr_dom = mort2yr_dom, 
                               recur2yr_dom = recur2yr_dom,
                               mort2yr_notDom = mort2yr_notDom, 
                               recur2yr_notDom = recur2yr_notDom, 
                               amp2yr_dom = amp2yr_dom,
                               amp2yr_notDom = amp2yr_notDom,
                               censorRate = censorRate)
  out <- bind_rows(out,
                   bind_cols(currParamDf, resultsOneRun) %>%
                     mutate(l = l, .before = "mort2yr_dom"))
}
# }
toc()

out$mort2yr_dom <- parameterGrid[command_args[5], "mort2yr1"]
out$recur2yr_dom <- parameterGrid[command_args[5], "woundRecur1"]
out$mort2yr_notDom <- parameterGrid[command_args[5], "mort2yr0"]
out$recur2yr_notDom <- parameterGrid[command_args[5], "woundRecur0"]
out$amp2yr_dom <- as.numeric(command_args[4])
out$amp2yr_notDom <- as.numeric(command_args[3])

outFileName <- paste0("3_aim1bPower_n", studySize, "_dom", dominantRegime, 
                      "_", round(command_args[3]*100), "_", round(command_args[4]*100), 
                      "_", command_args[5], ".csv")

if(str_detect(getwd(), "resubmission")){
  write_csv(x = out, file = paste0("./2_Aim1b/2_pipeline/", outFileName))
  } else{
  write_csv(x = out, file = paste0("../2_pipeline/", outFileName))
}





