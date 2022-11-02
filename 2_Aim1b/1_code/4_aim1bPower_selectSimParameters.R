# Power analysis for aim 1b
# Author: Nikki Freeman
# Last modified: 5 May 2022
# Dependencies (packages): tidyverse, sandwich
# Dependencies (scripts): 2_generateSampleAim1b.R

# Read the selected simulation parameters --------------------------------------
simParameters <- readxl::read_xlsx("./2_Aim1b/1_code/simulationParameters3.xlsx")

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
designProbs <- rep(0.5, 5)
L <- 1000
set.seed(20191022)

# Dataframe for output 
out <- data.frame()

# Simulation -------------------------------------------------------------------
tic()
for(i in 1:nrow(simParameters)){
  print(i)
  studySize <- simParameters$n[i]
  mort2yr_dom <- simParameters$mort1[i]
  recur2yr_dom <- simParameters$woundRecur1[i]
  mort2yr_notDom <- simParameters$mort0[i]
  recur2yr_notDom <- simParameters$woundRecur0[i]
  amp2yr_dom <- simParameters$amp1[i]
  amp2yr_notDom <- simParameters$amp0[i]
  dominantRegime <- simParameters$dominantRegime[i]
  censorRate <- simParameters$censorRate[i]
  currParamDf <- data.frame(mort2yr_dom = mort2yr_dom, 
                            recur2yr_dom = recur2yr_dom,
                            mort2yr_notDom = mort2yr_notDom, 
                            recur2yr_notDom = recur2yr_notDom,
                            amp2yr_dom = amp2yr_dom,
                            amp2yr_notDom = amp2yr_notDom,
                            censorRate = censorRate, 
                            dominantRegime= dominantRegime)
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
}
toc()

# out <- bind_cols(simParameters, out)

outFileName <- "4_aim1bPower_selectSimParameters3.csv"

if(str_detect(getwd(), "resubmission")){
  write_csv(x = out, file = paste0("./2_Aim1b/2_pipeline/", outFileName))
  } else{
  write_csv(x = out, file = paste0("../2_pipeline/", outFileName))
}





