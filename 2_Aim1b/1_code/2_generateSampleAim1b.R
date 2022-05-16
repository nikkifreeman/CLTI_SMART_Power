# Generate simulation data
# Author: Nikki Freeman
# Last modified: 4 May 2022
# Dependencies (packages): tidyverse

# Packages ---------------------------------------------------------------------
library(tidyverse)


# Generate patients for each trajectory (trajectory, or pathway, 1-11) ---------
# Draw response probs
generateObserved <- function(studySize, designProbs, responseProbRange,
                             dominantRegime, mort2yr_dom, 
                             recur2yr_dom, amp2yr_dom,
                             mort2yr_notDom, recur2yr_notDom, 
                             amp2yr_notDom, censorRate){
  # Draw the response probability
  p1 <- runif(1, min = responseProbRange$p1[1], max = responseProbRange$p1[2])
  p2 <- runif(1, min = responseProbRange$p2[1], max = responseProbRange$p2[2])
  p3 <- runif(1, min = responseProbRange$p3[1], max = responseProbRange$p3[2])
  p4 <- runif(1, min = responseProbRange$p4[1], max = responseProbRange$p4[2])
  p5 <- runif(1, min = responseProbRange$p5[1], max = responseProbRange$p5[2])
  
  # Generate the number in each trajectory -------------------------------------
  nA1_A <- round((studySize*designProbs[1])) # Number initially randomized to A
  nA1_B <- studySize - nA1_A # Number initially randomzed to B
  nA1_A_R <- rbinom(n = 1, size = nA1_A, prob = p1) # (1) Number of responders among those initially randomized to A
  nA1_A_NR <- nA1_A - nA1_A_R # Number of non-responders among those initially randomized to A
  nA1_A_A2_A <- round(nA1_A_NR * designProbs[2]) # Number randomized to A2 = A after failing A1 = A
  nA1_A_A2_B <- nA1_A_NR - nA1_A_A2_A # Number randomized to A2 = B after failing A1 = A
  nA1_A_A2_A_R <- rbinom(n = 1, size = nA1_A_A2_A, prob = p2) # (2) Number randomized to A2 = A after failing A1 = A that respond 
  nA1_A_A2_A_NR <- nA1_A_A2_A - nA1_A_A2_A_R # Number randomized to A2 = B after failing A1 = A that do not respond
  nA1_A_A2_A_B_R <- rbinom(n = 1, size = nA1_A_A2_A_NR, p = p3) # (3) Number A --> NR --> A --> B that respond
  nA1_A_A2_A_B_NR <- nA1_A_A2_A_NR - nA1_A_A2_A_B_R # Number A --> NR ---> A --> B that do not respond
  nA1_A_A2_A_B_A3_A <- round(nA1_A_A2_A_B_NR * designProbs[3]) # (4) Number A --> NR --> A --> NR --> B that are randomized to A3 = A
  nA1_A_A2_A_B_A3_D <- nA1_A_A2_A_B_NR - nA1_A_A2_A_B_A3_A # (5) Number A --> NR --> A --> NR --> B that are randomzied to A3 = D
  nA1_A_A2_B_R <- rbinom(n = 1, size = nA1_A_A2_B, prob = p4) # (6) Number A --> NR --> B that respond
  nA1_A_A2_B_NR <- nA1_A_A2_B - nA1_A_A2_B_R # Number A --> NR --> B that do not respond
  nA1_A_A2_B_A3_A <- round(nA1_A_A2_B_NR * designProbs[4]) # (7) Number A --> NR --> B --> NR that are randomized to A3 = A
  nA1_A_A2_B_A3_D <- nA1_A_A2_B_NR - nA1_A_A2_B_A3_A # (8) Number A --> NR --> B --> NR that are randomized to A3 = D
  nA1_B_R <- rbinom(n = 1, size = nA1_B, prob = p5) # (9) Number of responders among those initially randomized to B
  nA1_B_NR <- nA1_B - nA1_B_R # Number of non-responders among those initially radnomized to B
  nA1_B_A2_A <- round(nA1_B_NR * designProbs[5]) # (10) Number B --> NR that are randomized to A
  nA1_B_A2_D <- nA1_B_NR - nA1_B_A2_A # (11) Number B --> NR that are randomized to D
  
  # Data frames for the history of each trajectory -----------------------------
  trajectory1 <- data.frame(trajectory = rep(1, nA1_A_R),
                            # Treatment histories
                            A1 = "A", A2 = "none", A3 = "none",
                            # Regime consistency indicators
                            R1 = 1, R2 = 1, R3 = 1, R4 = 1, R5 = 0, R6 = 0)
  trajectory2 <- data.frame(trajectory = rep(2, nA1_A_A2_A_R ),
                            # Treatment histories
                            A1 = "A", A2 = "A", A3 = "none",
                            # Regime consistency indicators
                            R1 = 1, R2 = 1, R3 = 0, R4 = 0, R5 = 0, R6 = 0)
  trajectory3 <- data.frame(trajectory = rep(3, nA1_A_A2_A_B_R),
                            # Treatment histories
                            A1 = "A", A2 = "B", A3 = "none",
                            # Regime consistency indicators
                            R1 = 1, R2 = 1, R3 = 0, R4 = 0, R5 = 0, R6 = 0)
  trajectory4 <- data.frame(trajectory = rep(4, nA1_A_A2_A_B_A3_A),
                            # Treatment histories
                            A1 = "A", A2 = "A", A3 = "A",
                            # Regime consistency indicators
                            R1 = 1, R2 = 0, R3 = 0, R4 = 0, R5 = 0, R6 = 0)
  trajectory5 <- data.frame(trajectory = rep(5, nA1_A_A2_A_B_A3_D ),
                            # Treatment histories
                            A1 = "A", A2 = "A", A3 = "D",
                            # Regime consistency indicators
                            R1 = 0, R2 = 1, R3 = 0, R4 =0, R5 = 0, R6 = 0)
  trajectory6 <- data.frame(trajectory = rep(6, nA1_A_A2_B_R),
                            # Treatment histories
                            A1 = "A", A2 = "B", A3 = "none",
                            # Regime consistency indicators
                            R1 = 0, R2 = 0, R3 = 1, R4 = 1, R5 = 0, R6 = 0)
  trajectory7 <- data.frame(trajectory = rep(7, nA1_A_A2_B_A3_A),
                            # Treatment histories
                            A1 = "A", A2 = "B", A3 = "A",
                            # Regime consistency indicators
                            R1 = 0, R2 = 0, R3 = 1, R4 = 0, R5 = 0, R6 = 0)
  trajectory8 <- data.frame(trajectory = rep(8, nA1_A_A2_B_A3_D),
                            # Treatment histories
                            A1 = "A", A2 = "B", A3 = "D",
                            # Regime consistency indicators
                            R1 = 0, R2 = 0, R3 = 0, R4 = 1, R5 = 0, R6 = 0)
  trajectory9 <- data.frame(trajectory = rep(9, nA1_B_R),
                            # Treatment histories
                            A1 = "B", A2 = "none", A3 = "none",
                            # Regime consistency indicators
                            R1 = 0, R2 = 0, R3 = 0, R4 = 0, R5 = 1, R6 = 1)
  trajectory10 <- data.frame(trajectory = rep(10, nA1_B_A2_A),
                             # Treatment histories
                             A1 = "B", A2 = "A", A3 = "none",
                             # Regime consistency indicators
                             R1 = 0, R2 = 0, R3 = 0, R4 = 0, R5 = 1, R6 = 0)
  trajectory11 <- data.frame(trajectory = rep(11, nA1_B_A2_D),
                             # Treatment histories
                             A1 = "B", A2 = "D", A3 = "none",
                             # Regime consistency indicators
                             R1 = 0, R2 = 0, R3 = 0, R4 = 0, R5 = 0, R6 = 1)
  trajectories <- bind_rows(trajectory1, trajectory2, trajectory3,
                            trajectory4, trajectory5, trajectory6,
                            trajectory7, trajectory8, trajectory9, 
                            trajectory10, trajectory11)
  if(nrow(trajectories) != studySize){
    warning("Generated sample size does not match target study size")
  }
  
  # Generate outcomes for the trajectories consistent with the dominant regime ----
  trajectoriesConsistWithDomRegime <- trajectories %>%
    filter(!!sym(paste0("R", dominantRegime)) == 1)
  
  # Death is consistent for every trajectory, so calculate time of death for everyone
  timeToDeath_dom <- 365*2*(log(runif(n = nrow(trajectoriesConsistWithDomRegime)))/log(1-mort2yr_dom))
  trajectoriesConsistWithDomRegime <- trajectoriesConsistWithDomRegime %>%
    mutate(timeToDeath = timeToDeath_dom)
  
  # Amputation is consistent for every trajectory, so calculate time of amp for everyone
  timeToAmp_dom <- 365*2*(log(runif(n = nrow(trajectoriesConsistWithDomRegime)))/log(1-amp2yr_dom))
  trajectoriesConsistWithDomRegime <- trajectoriesConsistWithDomRegime %>%
    mutate(timeToAmp = timeToAmp_dom)
  
  # Wound recurrence is not consistent with trajectories 1 and 9
  trajectoriesConsistWithDomRegime_neverRecur <- trajectoriesConsistWithDomRegime %>%
    filter((trajectory %in% c(1, 9))) %>%
    mutate(timeToRecur = 10000, #These never recur so set recurrence far into the future
           timeToHealSecondWound = 0) %>% # These never recur so set healing time to 0
    rowwise() %>%
    mutate(timeToHealFirstWound = truncnorm::rtruncnorm(1, a = 0.5*30, mean = 2.5*30, sd = 2*30)) %>%
    ungroup() %>%
    select(trajectory, A1, A2, A3, R1, R2, R3, R4, R5, R6, timeToDeath, timeToAmp,
           timeToRecur, timeToHealFirstWound, timeToHealSecondWound) 
  
  # Deal with the trajectories that can have a recurrence
  trajectoriesConsistWithDomRegime_canRecur <- trajectoriesConsistWithDomRegime %>%
    # Trajectories other than 1 and 9 can have recurrences
    filter(!(trajectory %in% c(1, 9)))
  
  # Time to recurrence for those that can recur
  timeToRecur_dom <- 365*2*(log(runif(n = nrow(trajectoriesConsistWithDomRegime_canRecur)))/log(1-recur2yr_dom))
  
  # Get time to heal
  trajectoriesConsistWithDomRegime_canRecur <- trajectoriesConsistWithDomRegime_canRecur %>%
    mutate(timeToRecur = timeToRecur_dom) %>%
    rowwise() %>%
    mutate(timeToHealFirstWound = truncnorm::rtruncnorm(1, a = 0.5*30, mean = 2.5*30, sd = 2*30)) %>%
    mutate(timeToHealSecondWound = truncnorm::rtruncnorm(1, a = 0.5*30, mean = 2.5*30, sd = 2*30)) %>%
    ungroup() %>%
    select(trajectory, A1, A2, A3, R1, R2, R3, R4, R5, R6, timeToDeath, timeToAmp,
           timeToRecur, timeToHealFirstWound, timeToHealSecondWound)
  # Combine trajectories consistent with dominant regime
  trajectoriesConsistWithDomRegime <- bind_rows(trajectoriesConsistWithDomRegime_neverRecur, 
                                                trajectoriesConsistWithDomRegime_canRecur)
  
  # Generate outcomes for the trajectories not consistent with the dominant regime -----
  trajectoriesNOTConsistWithDomRegime <- trajectories %>%
    filter(!!sym(paste0("R", dominantRegime)) == 0)
  
  # Death is consistent for every trajectory, so calculate time of death for everyone
  timeToDeath_notDom <- 365*2*(log(runif(n = nrow(trajectoriesNOTConsistWithDomRegime)))/log(1-mort2yr_notDom))
  trajectoriesNOTConsistWithDomRegime <- trajectoriesNOTConsistWithDomRegime %>%
    mutate(timeToDeath = timeToDeath_notDom)
  
  # Amputation is consistent for every trajectory, so calculate time of amp for everyone
  timeToAmp_notDom <- 365*2*(log(runif(n = nrow(trajectoriesNOTConsistWithDomRegime)))/log(1-amp2yr_notDom))
  trajectoriesNOTConsistWithDomRegime <- trajectoriesNOTConsistWithDomRegime %>%
    mutate(timeToAmp = timeToAmp_notDom)
  
  # Wound recurrence is not consistent with trajectories 1 and 9
  trajectoriesNOTConsistWithDomRegime_neverRecur <- trajectoriesNOTConsistWithDomRegime %>%
    filter((trajectory %in% c(1, 9))) %>%
    mutate(timeToRecur = 10000, #These never recur so set recurrence far into the future
           timeToHealSecondWound = 0) %>% # These never recur so set healing time to 0
    rowwise() %>%
    mutate(timeToHealFirstWound = truncnorm::rtruncnorm(1, a = 0.5*30, mean = 2.5*30, sd = 2*30)) %>%
    ungroup() %>%
    select(trajectory, A1, A2, A3, R1, R2, R3, R4, R5, R6, timeToDeath, timeToAmp,
           timeToRecur, timeToHealFirstWound, timeToHealSecondWound) 
  
  # Deal with the trajectories that can have a recurrence
  trajectoriesNOTConsistWithDomRegime_canRecur <- trajectoriesNOTConsistWithDomRegime %>%
    # Trajectories other than 1 and 9 can have recurrences
    filter(!(trajectory %in% c(1, 9)))
  
  # Time to recurrence for those that can recur
  timeToRecur_notDom <- 365*2*(log(runif(n = nrow(trajectoriesNOTConsistWithDomRegime_canRecur)))/log(1-recur2yr_notDom))
  
  # Get time to heal
  trajectoriesNOTConsistWithDomRegime_canRecur <- trajectoriesNOTConsistWithDomRegime_canRecur %>%
    mutate(timeToRecur = timeToRecur_notDom) %>%
    rowwise() %>%
    mutate(timeToHealFirstWound = truncnorm::rtruncnorm(1, a = 0.5*30, mean = 2.5*30, sd = 2*30)) %>%
    mutate(timeToHealSecondWound = truncnorm::rtruncnorm(1, a = 0.5*30, mean = 2.5*30, sd = 2*30)) %>%
    ungroup() %>%
    select(trajectory, A1, A2, A3, R1, R2, R3, R4, R5, R6, timeToDeath, timeToAmp,
           timeToRecur, timeToHealFirstWound, timeToHealSecondWound)
  
  # Combine trajectories consistent with dominant regime
  trajectoriesNOTConsistWithDomRegime <- bind_rows(trajectoriesNOTConsistWithDomRegime_neverRecur, 
                                                trajectoriesNOTConsistWithDomRegime_canRecur)
  
  # Combine get wound-free alive days and add censoring ------------------------
  trajectoriesWithTTE <- bind_rows(trajectoriesConsistWithDomRegime, trajectoriesNOTConsistWithDomRegime) %>%
    mutate(firstWoundHealingDay = timeToHealFirstWound) %>%
    rowwise() %>%
    mutate(timeToRecur = max(timeToRecur, firstWoundHealingDay + 1)) %>%
    ungroup() %>%
    mutate(secondWoundHealingDay = timeToRecur + timeToHealSecondWound) %>%
    mutate(endOfStudy = 365*2) %>%
    mutate(id = 1:studySize, .before = "trajectory")  %>%
    mutate(ltfu = if_else(id %in% sample(1:studySize, size = round(studySize*censorRate), replace = FALSE), 1, 0))
  
  # for convenience, make event1:event5 variables
  sequences <- trajectoriesWithTTE %>%
    pivot_longer(cols = c(timeToDeath, timeToAmp, timeToRecur,
                          firstWoundHealingDay, secondWoundHealingDay, endOfStudy),
                 names_to = "events", values_to = "day") %>%
    group_by(id) %>% arrange(id, day) %>% select(id, events, day) %>% ungroup() %>%
    group_by(id) %>%
    summarise(sequence = paste(events, collapse = "-")) %>%
    mutate(event1 = str_extract(sequence, "^\\w*(?=-)")) %>%
    mutate(sequence_trimmed = str_remove(sequence, "^\\w*-")) %>%
    mutate(event2 = str_extract(sequence_trimmed, "^\\w*(?=-)")) %>%
    mutate(sequence_trimmed = str_remove(sequence_trimmed, "^\\w*-")) %>%
    mutate(event3 = str_extract(sequence_trimmed, "\\w*(?=-)")) %>%
    mutate(sequence_trimmed = str_remove(sequence_trimmed, "^\\w*-")) %>%
    mutate(event4 = str_extract(sequence_trimmed, "\\w*(?=-)")) %>%
    mutate(sequence_trimmed = str_remove(sequence_trimmed, "^\\w*-")) %>%
    mutate(event5 = str_extract(sequence_trimmed, "\\w*(?=-)")) %>%
    mutate(sequence_trimmed = str_remove(sequence_trimmed, "\\w*-")) %>%
    mutate(event6 = str_extract(sequence_trimmed, "^\\w*")) %>%
    ungroup() %>%
    select(-sequence_trimmed)
  
  
  allDf <- left_join(trajectoriesWithTTE, sequences, by = "id")
  
  allDf %>%
    # Characterize the outcomes
    # The first event is one of the following: heal the first wound, die, amp, or end of study
    mutate(outcome_qual = if_else(event1 == "timeToDeath",
                                  "died before healing first wound and before study end",
                                  "?")) %>%
    mutate(outcome_qual = if_else(event1 == "eos",
                                  "did not heal first wound, alive at study end",
                                  outcome_qual)) %>%
    mutate(outcome_qual = if_else(event1 == "timeToAmp" & timeToDeath <= endOfStudy,
                                  "amp before healing first wound then died before study end",
                                  outcome_qual)) %>%
    mutate(outcome_qual = if_else(event1 == "timeToAmp" & timeToDeath > endOfStudy,
                                  "amp before healing first wound, alive at study end",
                                  outcome_qual)) %>%
    # Of those that healed first wound, the second event is one of four events: die, amp, recur, or end of study
    mutate(outcome_qual = if_else(event1 == "firstWoundHealingDay" &
                                    event2 == "timeToDeath",
                                  "healed first wound, died before recur and before study end",
                                  outcome_qual)) %>%
    mutate(outcome_qual = if_else(event1 == "firstWoundHealingDay" &
                                    event2 == "endOfStudy",
                                  "healed first wound, no recur, alive with limb at study end",
                                  outcome_qual)) %>%
    mutate(outcome_qual = if_else(event1 == "firstWoundHealingDay" &
                                    event2 == "timeToAmp" &
                                    timeToDeath <= endOfStudy,
                                  "healed first wound, amp, died before study end",
                                  outcome_qual)) %>%
    mutate(outcome_qual = if_else(event1 == "firstWoundHealingDay" &
                                    event2 == "timeToAmp" &
                                    timeToDeath > endOfStudy,
                                  "healed first wound, amp, alive at study end",
                                  outcome_qual)) %>%
    # Of those that recurred, the third event can only be one of: die, amp, heal, end of study
    mutate(outcome_qual = if_else(event1 == "firstWoundHealingDay" &
                                    event2 == "timeToRecur" &
                                    event3 == "timeToDeath",
                                  "healed first, recurred, died before healing second",
                                  outcome_qual)) %>%
    mutate(outcome_qual = if_else(event1 == "firstWoundHealingDay" &
                                    event2 == "timeToRecur" &
                                    event3 == "endOfStudy",
                                  "healed, recurred, end of study before healing 2nd wound",
                                  outcome_qual)) %>%
    mutate(outcome_qual = if_else(event1 == "firstWoundHealingDay" &
                                    event2 == "timeToRecur" &
                                    event3 == "timeToAmp" &
                                    timeToDeath <= endOfStudy,
                                  "healed, recurred, amp before healing 2nd, dead at study end",
                                  outcome_qual)) %>%
    mutate(outcome_qual = if_else(event1 == "firstWoundHealingDay" &
                                    event2 == "timeToRecur" &
                                    event3 == "timeToAmp" &
                                    timeToDeath > endOfStudy,
                                  "healed, recurred, amp before healing 2nd, alive at study end",
                                  outcome_qual)) %>%
    # Of those the healed their recurrent wound, event 4 can only be one of: die, end of study, amp
    mutate(outcome_qual = if_else(event1 == "firstWoundHealingDay" &
                                    event2 == "timeToRecur" &
                                    event3 == "secondWoundHealingDay" &
                                    event4 == "timeToDeath",
                                  "healed, recurred, healed, died",
                                  outcome_qual)) %>%
    mutate(outcome_qual = if_else(event1 == "firstWoundHealingDay" &
                                    event2 == "timeToRecur" &
                                    event3 == "secondWoundHealingDay" &
                                    event4 == "endOfStudy",
                                  "healed, recurred, healed, alive with limb at study end",
                                  outcome_qual)) %>%
    mutate(outcome_qual = if_else(event1 == "firstWoundHealingDay" &
                                    event2 == "timeToRecur" &
                                    event3 == "secondWoundHealingDay" &
                                    event4 == "timeToAmp" &
                                    timeToDeath <= endOfStudy,
                                  "healed, recurred, healed, amp, died",
                                  outcome_qual)) %>%
    mutate(outcome_qual = if_else(event1 == "firstWoundHealingDay" &
                                    event2 == "timeToRecur" &
                                    event3 == "secondWoundHealingDay" &
                                    event4 == "timeToAmp" &
                                    timeToDeath > endOfStudy,
                                  "healed, recurred, healed, amp, alive at study end",
                                  outcome_qual)) %>%
    # Calculate wound free alive days
    mutate(woundFreeAliveDays = if_else(outcome_qual == "died before healing first wound and before study end",
                                        0, -1000)) %>%
    mutate(woundFreeAliveDays = if_else(outcome_qual == "did not heal first wound, alive at study end",
                                        0, woundFreeAliveDays)) %>%
    mutate(woundFreeAliveDays = if_else(outcome_qual == "amp before healing first wound then died before study end",
                                        timeToDeath - timeToAmp, woundFreeAliveDays)) %>%
    mutate(woundFreeAliveDays = if_else(outcome_qual == "amp before healing first wound, alive at study end",
                                        endOfStudy - timeToAmp, woundFreeAliveDays)) %>%
    mutate(woundFreeAliveDays = if_else(outcome_qual == "healed first wound, died before recur and before study end",
                                        timeToDeath - firstWoundHealingDay, woundFreeAliveDays)) %>%
    mutate(woundFreeAliveDays = if_else(outcome_qual == "healed first wound, no recur, alive with limb at study end",
                                        endOfStudy - firstWoundHealingDay, woundFreeAliveDays)) %>%
    mutate(woundFreeAliveDays = if_else(outcome_qual == "healed first wound, amp, died before study end",
                                        timeToDeath - firstWoundHealingDay, woundFreeAliveDays)) %>%
    mutate(woundFreeAliveDays = if_else(outcome_qual == "healed first wound, amp, alive at study end",
                                        endOfStudy - firstWoundHealingDay, woundFreeAliveDays)) %>%
    mutate(woundFreeAliveDays = if_else(outcome_qual == "healed first, recurred, died before healing second",
                                        timeToRecur - firstWoundHealingDay, woundFreeAliveDays)) %>%
    mutate(woundFreeAliveDays = if_else(outcome_qual == "healed, recurred, end of study before healing 2nd wound",
                                        timeToRecur - firstWoundHealingDay, woundFreeAliveDays)) %>%
    mutate(woundFreeAliveDays = if_else(outcome_qual == "healed, recurred, amp before healing 2nd, dead at study end",
                                        (timeToRecur - firstWoundHealingDay) + (timeToDeath - timeToAmp), woundFreeAliveDays)) %>%
    mutate(woundFreeAliveDays = if_else(outcome_qual == "healed, recurred, amp before healing 2nd, alive at study end",
                                        (timeToRecur - firstWoundHealingDay) + (endOfStudy - timeToAmp), woundFreeAliveDays)) %>%
    mutate(woundFreeAliveDays = if_else(outcome_qual == "healed, recurred, healed, died",
                                        (timeToRecur - firstWoundHealingDay) + (timeToDeath - secondWoundHealingDay), woundFreeAliveDays)) %>%
    mutate(woundFreeAliveDays = if_else(outcome_qual == "healed, recurred, healed, alive with limb at study end",
                                        (timeToRecur - firstWoundHealingDay) + (endOfStudy - secondWoundHealingDay), woundFreeAliveDays)) %>%
    mutate(woundFreeAliveDays = if_else(outcome_qual == "healed, recurred, healed, amp, died",
                                        (timeToRecur - firstWoundHealingDay) + (timeToDeath - secondWoundHealingDay), woundFreeAliveDays)) %>%
    mutate(woundFreeAliveDays = if_else(outcome_qual == "healed, recurred, healed, amp, alive at study end",
                                        (timeToRecur - firstWoundHealingDay) + (endOfStudy - secondWoundHealingDay), woundFreeAliveDays))
    
    

}



