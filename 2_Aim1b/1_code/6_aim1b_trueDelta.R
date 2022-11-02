# Compute the true deltas
# Author: Nikki Freeman
# Last modified: 19 June 2022
# Dependencies (packages): tidyverse, tictoc
# Dependenices: 2_generateSampleAim1b.R, simulationParameters.xlsx

# Packages
library(tidyverse)
library(tictoc)

# Load the required functions and the simulation parameters
if(str_detect(getwd(), "resub")){
  source("./2_Aim1b/1_code/2_generateSampleAim1b.R")
  simParameters <- readxl::read_xlsx("./2_Aim1b/1_code/simulationParameters3.xlsx")
} else{
  source("2_generateSampleAim1b.R")
  simParameters <- readxl::read_xlsx("simulationParameters3.xlsx")
}

# Placeholder for results
out <- data.frame()

# Simulation parameters 
set.seed(5789)
responseProbRange <- list(p1 = c(0.25, 0.35),
                          p2 = c(0.25, 0.35),
                          p3 = c(0.3, 0.6),
                          p4 = c(0.3, 0.6),
                          p5 = c(0.3, 0.6))
tic()
for(i in 1:nrow(simParameters)){
  temp <- generateObserved(studySize = 10000, 
                           designProbs = rep(0.5, 5), 
                           responseProbRange = responseProbRange, 
                           dominantRegime = simParameters$dominantRegime[i], 
                           mort2yr_dom = simParameters$mort1[i], 
                           recur2yr_dom = simParameters$woundRecur1[i], 
                           amp2yr_dom = simParameters$amp1[i], 
                           mort2yr_notDom = simParameters$mort0[i], 
                           recur2yr_notDom = simParameters$woundRecur0[i], 
                           amp2yr_notDom = simParameters$amp0[i], 
                           censorRate = simParameters$censorRate[i])
  mu_dom <- temp %>% filter(!!sym(paste0("R", simParameters$dominantRegime[i])) == 1) %>%
    summarise(mean1 = mean(woundFreeAliveDays), sd1 = sd(woundFreeAliveDays))
  mu_notDom <- temp %>% filter(!!sym(paste0("R", simParameters$dominantRegime[i])) != 1) %>%
    summarise(mean0 = mean(woundFreeAliveDays), sd0 = sd(woundFreeAliveDays))
  
  toOutput_df <- data.frame(dominantRegime = simParameters$dominantRegime[i], 
                            mort2yr_dom = simParameters$mort1[i], 
                            recur2yr_dom = simParameters$woundRecur1[i], 
                            amp2yr_dom = simParameters$amp1[i], 
                            mort2yr_notDom = simParameters$mort0[i], 
                            recur2yr_notDom = simParameters$woundRecur0[i], 
                            amp2yr_notDom = simParameters$amp0[i], 
                            censorRate = simParameters$censorRate[i],
                            mu_dom = mu_dom$mean1,
                            sd_dom = mu_dom$sd1,
                            mu_notdom = mu_notDom$mean0,
                            sd_notDom = mu_notDom$sd0)
  
  out <- bind_rows(out, toOutput_df)
}
toc()

# save the results
if(str_detect(getwd(), "resub")){
  write_csv(x = out, file = "./2_Aim1b/2_pipeline/6_aim1b_trueDelta3.csv")
} else{
  write_csv(x = out, file = "../2_pipeline/6_aim1b_trueDelta3.csv")
}
