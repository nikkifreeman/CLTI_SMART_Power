# 0_universalParameters.R
# Author: Nikki Freeman
# Last modified: 9 May 2022

# Preamble ---------------------------------------------------------------------
# This file contains the universal simulation parameters used for all of the 
# power simulations
library(gtools)

# Simulation parameters 
# N <- c(550, 600, 650, 700, 750) # Study size
C_L <- 0.15 # LTFU rate
mort2yr <- seq(from = 0.15, to = 0.30, by = 0.05)
woundRecur <- seq(from = 0.25, to = 0.45, by = 0.05)
amp2yr <- c(0.20, 0.25)
L <- 1000 # Number of simulations for power calc
responseProbRange <- list(p1 = c(0.25, 0.35),
                          p2 = c(0.25, 0.35),
                          p3 = c(0.3, 0.6),
                          p4 = c(0.3, 0.6),
                          p5 = c(0.3, 0.6))

# Create a grid of 2-year mort and 2-year wound recurrence values to evaluate over
mortalityCombinations <- gtools::permutations(n = length(mort2yr), r = 2, v = mort2yr, repeats.allowed = TRUE)
woundRecurCombinations <- gtools::permutations(n = length(woundRecur), r = 2, v = woundRecur, repeats.allowed = TRUE)
parameterGrid <- data.frame()
for(i in 1:nrow(mortalityCombinations)){
  newParams <- data.frame(mort2yr0 = mortalityCombinations[i, 1], mort2yr1 = mortalityCombinations[i, 2]) %>%
    cbind(woundRecurCombinations)
  parameterGrid <- bind_rows(parameterGrid, newParams)
}
names(parameterGrid) <- c("mort2yr0", "mort2yr1", "woundRecur0", "woundRecur1")
parameterGrid <- parameterGrid %>% mutate(C_L = C_L)

# Create a grid of 2-year mort, 2-year wound recurrence, and 2-year amputation values
# to evaluate over
# parameterGrid_withAmp <- expand.grid(mort2yr0 = mort2yr, mort2yr1 = mort2yr,
#                                      woundRecur0 = woundRecur, woundRecur1 = woundRecur,
#                                      amp2yr0 = amp2yr, amp2yr1 = amp2yr) %>%
#   filter(!(mort2yr1 == mort2yr0 & woundRecur1 == woundRecur0 & amp2yr1 == amp2yr0))

parameterGrid_withAmp <- parameterGrid %>%
  mutate(amp0 = 0.2, amp1 = 0.2, .after = woundRecur1) %>%
  bind_rows(parameterGrid %>%
              mutate(amp0 = 0.2, amp1 = 0.25, .after = woundRecur1)) %>%
  bind_rows(parameterGrid %>%
              mutate(amp0 = 0.25, amp1 = 0.2, .after = woundRecur1)) %>%
  bind_rows(parameterGrid %>% 
              mutate(amp0 = 0.25, amp1 = 0.25, .after = woundRecur1))


