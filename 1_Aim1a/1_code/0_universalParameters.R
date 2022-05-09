# 0_universalParameters.R
# Author: Nikki Freeman
# Last modified: 25 April 2022

# Preamble ---------------------------------------------------------------------
# This file contains the universal simulation parameters used for all of the 
# power simulations

# Simulation parameters 
N <- c(550, 600, 650, 700, 750) # Study size
C_L <- 0.15 # LTFU rate
mort2yr <- seq(from = 0.15, to = 0.45, by = 0.05)
woundRecur <- seq(from = 0.25, to = 0.55, by = 0.05)
L <- 1000 # Number of simulations for power calc
responseProbRange <- list(p1 = c(0.25, 0.35),
                          p2 = c(0.25, 0.35),
                          p3 = c(0.3, 0.6),
                          p4 = c(0.3, 0.6),
                          p5 = c(0.3, 0.6))

# Create a grid of 2-year mort and 2-year wound recurrence values to evaluate over
parameterGrid <- expand.grid(mort2yr0 = mort2yr, mort2yr1 = mort2yr, woundRecur0 = woundRecur, woundRecur1 = woundRecur, C_L = C_L) %>%
  # Remove the rows where mortality and wound recurrence are the same in both arms
  filter(!(mort2yr1 == mort2yr0 & woundRecur1 == woundRecur0)) 