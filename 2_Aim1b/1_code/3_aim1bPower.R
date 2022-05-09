# Power analysis for aim 1b
# Author: Nikki Freeman
# Last modified: 5 May 2022
# Dependencies (packages): tidyverse, sandwich
# Dependencies (scripts): 2_generateSampleAim1b.R

# Packages ---------------------------------------------------------------------
library(tidyverse)
library(sandwich)

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
studySize <- 650
designProbs <- rep(0.5, 5)
dominantRegime <- 1

# Dataframe for output 
out <- data.frame()

# Simulation -------------------------------------------------------------------
for(i in 1:nrow(parameterGrid)){
  mort2yr_dom <- parameterGrid[i, "mort2yr1"]
  recur2yr_dom <- parameterGrid[i, "woundRecur1"]
  mort2yr_notDom <- parameterGrid[i, "mort2yr0"]
  recur2yr_notDom <- parameterGrid[i, "woundRecur0"]
  censorRate <- C_L
  currParamDf <- data.frame(mort2yr_dom = mort2yr_dom, 
                            recur2yr_dom = recur2yr_dom,
                            mort2yr_notDom = mort2yr_notDom, 
                            recur2yr_notDom = recur2yr_notDom,
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
                                 censorRate = censorRate)
    out <- bind_rows(out,
              bind_cols(currParamDf, resultsOneRun) %>%
                mutate(l = l, .before = "mort2yr_dom"))
  }
}

outFileName <- paste0("3_aim1bPower_n", studySize, "_dom", dominantRegime, ".csv")

if(str_detect(getwd(), "resubmission")){
  write_csv(x = out, file = paste0("./2_Aim1b/2_pipeline/", outFileName))
} else{
  write_csv(x = out, file = paste0("../2_pipeline/", outFileName))
}




