# 1_generateSimulationData.R
# Author: Nikki Freeman
# Last modified: 16 May 2022

# Packages ---------------------------------------------------------------------
library(tidyverse)

# Generate trajectories and outcomes -------------------------------------------
generateObserved <- function(studySize, designProbs, responseProbRange,
                             dominantRegime, mort2yr_dom, 
                             recur2yr_dom, amp2yr_dom,
                             mort2yr_notDom, recur2yr_notDom, 
                             amp2yr_notDom, censorRate){
  # Draw the response probabilities (these will be the probabilities conditional on X_n = 1 and X_n = 0)
  p1 <- runif(2, min = responseProbRange$p1[1], max = responseProbRange$p1[2])
  p2 <- runif(2, min = responseProbRange$p2[1], max = responseProbRange$p2[2])
  p3 <- runif(1, min = responseProbRange$p3[1], max = responseProbRange$p3[2])
  p4 <- runif(2, min = responseProbRange$p4[1], max = responseProbRange$p4[2])
  p5 <- runif(2, min = responseProbRange$p5[1], max = responseProbRange$p5[2])
  
  # Generate trajectories  -----------------------------------------------------
  # Generate the tailoring variables
  X1 <- rbinom(n = studySize, size = 1, prob = 0.4)
  X2 <- rbinom(n = studySize, size = 1, prob = 0.5)
  X3 <- rbinom(n = studySize, size = 1, prob = 0.6)
  X4 <- rbinom(n = studySize, size = 1, prob = 0.45)
  X5 <- rbinom(n = studySize, size = 1, prob = 0.55)
  # Generate the first treatment assignment
  A1 <- (1:studySize %in% sample(1:studySize, size = studySize*designProbs[1], replace = FALSE))*1.0 # Generate the first stage treatment assignments
  A1[A1 == 1] <- "A"
  A1[A1 == "0"] <- "B"
  # Generate the responses
  temp <- data.frame(id = seq(from = 1, to = studySize, by = 1), X1, X2, X3, X4, X5, A1) %>%
    mutate(R1 = -1)
  temp$R1[temp$A1 == 'A' & temp$X1 == 1] <- rbinom(n = nrow(temp[temp$A1 == 'A' & temp$X1 == 1,]), size = 1, max(p1))
  temp$R1[temp$A1 == "A" & temp$X1 == 0] <- rbinom(n = nrow(temp[temp$A1 == "A" & temp$X1 == 0,]), size = 1, min(p1))
  temp$R1[temp$A1 == "B" & temp$X1 == 0] <- rbinom(n = nrow(temp[temp$A1 == "B" & temp$X1 == 0,]), size = 1, max(p5))
  temp$R1[temp$A1 == "B" & temp$X1 == 1] <- rbinom(n = nrow(temp[temp$A1 == "B" & temp$X1 == 1,]), size = 1, min(p5))
  # Generate the 2nd treatment assignment
  temp <- temp %>% mutate(A2 = "?")
  temp$A2[temp$A1 == "B" & temp$R1 == 0] <- temp$id[temp$A1 == "B" & temp$R1 == 0] %in% sample(temp$id[temp$A1 == "B" & temp$R1 == 0], 
                                                     size = round(length(temp$id[temp$A1 == "B" & temp$R1 == 0])*designProbs[5]))
  temp$A2[temp$A2 == "TRUE"] <- "A"
  temp$A2[temp$A2 == "FALSE"] <- "D"
  temp$A2[temp$A1 == "A" & temp$R1 == 0] <- temp$id[temp$A1 == "A" & temp$R1 == 0] %in% sample(temp$id[temp$A1 == "A" & temp$R1 == 0],
                                                    size = round(length(temp$id[temp$A1 == 'A' & temp$R1 == 0])*designProbs[2]))
  temp$A2[temp$A2 == "TRUE"] <- "A"
  temp$A2[temp$A2 == "FALSE"] <- "B"
  temp$A2[temp$A2 == "?"] <- "first stage responder"
  # Generate the responses after the 2nd randomization
  temp$R2 <- -1
  temp$R2[temp$A1 == "A" & temp$A2 == "A" & temp$X2 == 1] <- rbinom(n = length(temp$R2[temp$A1 == "A" & temp$A2 == "A" & temp$X2 == 1]),
                                                                    size = 1, max(p2) + (1-p2)*p3)
  temp$R2[temp$A1 == "A" & temp$A2 == "A" & temp$X2 == 0] <- rbinom(n = length(temp$R2[temp$A1 == "A" & temp$A2 == "A" & temp$X2 == 0]),
                                                                    size = 1, min(p2) + (1-p2)*p3)
  temp$R2[temp$A1 == "A" & temp$A2 == "B" & temp$X2 == 0] <- rbinom(n = length(temp$R2[temp$A1 == "A" & temp$A2 == "B" & temp$X2 == 0]),
                                                                    size = 1, max(p4))
  temp$R2[temp$A1 == "A" & temp$A2 == "B" & temp$X2 == 1] <- rbinom(n = length(temp$R2[temp$A1 == "A" & temp$A2 == "B" & temp$X2 == 1]),
                                                                    size = 1, min(p4))
  # Generate the 3rd treatment assignment
  temp$A3 <- "?"
  temp$A3[temp$A2 == "A" & temp$R2 == 0] <- temp$id[temp$A2 == "A" & temp$R2 == 0] %in% sample(temp$id[temp$A2 == "A" & temp$R2 == 0],
                                                                                               size = round(length(temp$id[temp$A2 == "A" & temp$R2 == 0])*designProbs[3]))
  temp$A3[temp$A3 == "TRUE"] <- "A"
  temp$A3[temp$A3 == "FALSE"] <- "D"
  temp$A3[temp$A2 == "B" & temp$R2 == 0] <- temp$id[temp$A2 == "B" & temp$R2 == 0] %in% sample(temp$id[temp$A2 == "B" & temp$R2 == 0],
                                                                                               size = round(length(temp$id[temp$A2 == "B" & temp$R2 == 0])*designProbs[4]))
  temp$A3[temp$A3 == "TRUE"] <- "A"
  temp$A3[temp$A3 == "FALSE"] <- "D"
  temp$A3[temp$A3 == "?" & temp$R2 == 1] <- "second stage responder"
  temp$A3[temp$A1 == "B"] <- "not eligible for 3rd randomization (A1 == B)"
  temp$A3[temp$A2 == "first stage responder"] <- "first stage responder"
  
  # add the trajectory labels for convenience
  temp$trajectory <- "?"
  temp$trajectory[temp$A1 == "A" & temp$R1 == 1] <- 1
  temp$trajectory[temp$A1 == "A" & temp$R1 == 0 & temp$A2 == "A" & temp$R2 == 1] <- 23
  temp$trajectory[temp$A1 == "A" & temp$R1 == 0 & temp$A2 == "A" & temp$R2 == 0 & temp$A3 == "A"] <- 4
  temp$trajectory[temp$A1 == "A" & temp$R1 == 0 & temp$A2 == "A" & temp$R2 == 0 & temp$A3 == "D"] <- 5
  temp$trajectory[temp$A1 == "A" & temp$R1 == 0 & temp$A2 == "B" & temp$R2 == 1] <- 6
  temp$trajectory[temp$A1 == "A" & temp$R1 == 0 & temp$A2 == "B" & temp$R2 == 0 & temp$A3 == "A"] <- 7
  temp$trajectory[temp$A1 == "A" & temp$R1 == 0 & temp$A2 == "B" & temp$R2 == 0 & temp$A3 == "D"] <- 8
  temp$trajectory[temp$A1 == "B" & temp$R1 == 1] <- 9
  temp$trajectory[temp$A1 == "B" & temp$R1 == 0 & temp$A2 == "A"] <- 10
  temp$trajectory[temp$A1 == "B" & temp$R1 == 0 & temp$A2 == "D"] <- 11
  
  # Generate outcomes for trajectories  ----------------------------------------
  # S(t) = \exp{-\lambda t}
  # -1 * log(S(t))/t = \lambda
  # T = -log(u)/\lambda*exp(beta^T X)
  # #
  # Outcomes for those in trajectory 1 (function of X1)
  trajectory1 <- temp %>% 
    filter(trajectory == 1) %>%
    # Calculate time to death
    rowwise() %>%
    mutate(timeToDeath = if_else(X1 == 1, 365*2*(log(runif(n = 1))/log(1-mort2yr_dom)),
                                 365*2*(log(runif(n = 1))/log(1-mort2yr_notDom)))) %>%
    # Calculate time to amputation
    mutate(timeToAmp = if_else(X1 == 1, 365*2*(log(runif(n = 1))/log(1-amp2yr_dom)),
                               365*2*(log(runif(n = 1))/log(1-amp2yr_notDom)))) %>%
    # Calculate time to heal first wound
    mutate(timeToHealFirstWound = truncnorm::rtruncnorm(1, a = 0.5*30, mean = 2.5*30, sd = 2*30)) %>%
    ungroup() %>%
    mutate(timeToRecur = 10000) %>%
    mutate(timeToHealSecondWound = 0)
  
  # Outcomes for those in trajectory 2 or 3 (function of X1 and X2)
  trajectory23 <- temp %>% 
    filter(trajectory == 23) %>%
    # Calculate time to death 
    rowwise() %>%
    mutate(timeToDeath = if_else(X2 == 1, 365*2*(log(runif(n = 1))/log(1-mort2yr_dom)),
                                 365*2*(log(runif(n = 1))/log(1-mort2yr_notDom)))) %>%
    # Calculate the time to amputation
    mutate(timeToAmp = if_else(X2 == 1, 365*2*(log(runif(n = 1))/log(1-amp2yr_dom)),
                               365*2*(log(runif(n = 1))/log(1-amp2yr_notDom)))) %>%
    # Calculate time to heal first wound
    mutate(timeToHealFirstWound = truncnorm::rtruncnorm(1, a = 0.5*30, mean = 2.5*30, sd = 2*30)) %>%
    # Calculate time to recurrence
    mutate(timeToRecur = if_else(X2 == 1, 365*2*(log(runif(n = 1))/log(1-recur2yr_dom)),
                                 365*2*(log(runif(n = 1))/log(1-recur2yr_notDom)))) %>%
    # Calculate time to heal second wound
    mutate(timeToHealSecondWound = truncnorm::rtruncnorm(1, a = 0.5*30, mean = 2.5*30, sd = 2*30)) %>%
    ungroup() 

  
  # Outcomes for those in trajectory 4 (function of X1, X2, and X3)
  trajectory4 <- temp %>%
    filter(trajectory == 4) %>%
    # Calculate time to death 
    rowwise() %>%
    mutate(timeToDeath = if_else(X3 == 1, 365*2*(log(runif(n = 1))/log(1-mort2yr_dom)),
                                 365*2*(log(runif(n = 1))/log(1-mort2yr_notDom)))) %>%
    # Calculate the time to amputation
    mutate(timeToAmp = if_else(X3 == 1, 365*2*(log(runif(n = 1))/log(1-amp2yr_dom)),
                               365*2*(log(runif(n = 1))/log(1-amp2yr_notDom)))) %>%
    # Calculate time to heal first wound
    mutate(timeToHealFirstWound = truncnorm::rtruncnorm(1, a = 0.5*30, mean = 2.5*30, sd = 2*30)) %>%
    # Calculate time to recurrence
    mutate(timeToRecur = if_else(X3 == 1, 365*2*(log(runif(n = 1))/log(1-recur2yr_dom)),
                                 365*2*(log(runif(n = 1))/log(1-recur2yr_notDom)))) %>%
    # Calculate time to heal second wound
    mutate(timeToHealSecondWound = truncnorm::rtruncnorm(1, a = 0.5*30, mean = 2.5*30, sd = 2*30)) %>%
    ungroup() 
  
  # Outcomes for those in trajectory 5 (function of X1, X2, and X3)
  trajectory5 <- temp %>%
    filter(trajectory == 5) %>%
    # Calculate time to death 
    rowwise() %>%
    mutate(timeToDeath = if_else(X3 == 0, 365*2*(log(runif(n = 1))/log(1-mort2yr_dom)),
                                 365*2*(log(runif(n = 1))/log(1-mort2yr_notDom)))) %>%
    # Calculate the time to amputation
    mutate(timeToAmp = if_else(X3 == 0, 365*2*(log(runif(n = 1))/log(1-amp2yr_dom)),
                               365*2*(log(runif(n = 1))/log(1-amp2yr_notDom)))) %>%
    # Calculate time to heal first wound
    mutate(timeToHealFirstWound = truncnorm::rtruncnorm(1, a = 0.5*30, mean = 2.5*30, sd = 2*30)) %>%
    # Calculate time to recurrence
    mutate(timeToRecur = if_else(X3 == 0, 365*2*(log(runif(n = 1))/log(1-recur2yr_dom)),
                                 365*2*(log(runif(n = 1))/log(1-recur2yr_notDom)))) %>%
    # Calculate time to heal second wound
    mutate(timeToHealSecondWound = truncnorm::rtruncnorm(1, a = 0.5*30, mean = 2.5*30, sd = 2*30)) %>%
    ungroup() 
  
  # Outcomes for those in trajectory 6 (function of X1 and X2)
  trajectory6 <- temp %>%
    filter(trajectory == 6) %>%
    # Calculate time to death 
    rowwise() %>%
    mutate(timeToDeath = if_else(X2 == 0, 365*2*(log(runif(n = 1))/log(1-mort2yr_dom)),
                                 365*2*(log(runif(n = 1))/log(1-mort2yr_notDom)))) %>%
    # Calculate the time to amputation
    mutate(timeToAmp = if_else(X2 == 0, 365*2*(log(runif(n = 1))/log(1-amp2yr_dom)),
                               365*2*(log(runif(n = 1))/log(1-amp2yr_notDom)))) %>%
    # Calculate time to heal first wound
    mutate(timeToHealFirstWound = truncnorm::rtruncnorm(1, a = 0.5*30, mean = 2.5*30, sd = 2*30)) %>%
    # Calculate time to recurrence
    mutate(timeToRecur = if_else(X2 == 0, 365*2*(log(runif(n = 1))/log(1-recur2yr_dom)),
                                 365*2*(log(runif(n = 1))/log(1-recur2yr_notDom)))) %>%
    # Calculate time to heal second wound
    mutate(timeToHealSecondWound = truncnorm::rtruncnorm(1, a = 0.5*30, mean = 2.5*30, sd = 2*30)) %>%
    ungroup() 
  
  # Outcomes for those in trajectory 7 (function of X1, X2, X4)
  trajectory7 <- temp %>%
    filter(trajectory == 7) %>%
    # Calculate time to death 
    rowwise() %>%
    mutate(timeToDeath = if_else(X4 == 1, 365*2*(log(runif(n = 1))/log(1-mort2yr_dom)),
                                 365*2*(log(runif(n = 1))/log(1-mort2yr_notDom)))) %>%
    # Calculate the time to amputation
    mutate(timeToAmp = if_else(X4 == 1, 365*2*(log(runif(n = 1))/log(1-amp2yr_dom)),
                               365*2*(log(runif(n = 1))/log(1-amp2yr_notDom)))) %>%
    # Calculate time to heal first wound
    mutate(timeToHealFirstWound = truncnorm::rtruncnorm(1, a = 0.5*30, mean = 2.5*30, sd = 2*30)) %>%
    # Calculate time to recurrence
    mutate(timeToRecur = if_else(X4 == 1, 365*2*(log(runif(n = 1))/log(1-recur2yr_dom)),
                                 365*2*(log(runif(n = 1))/log(1-recur2yr_notDom)))) %>%
    # Calculate time to heal second wound
    mutate(timeToHealSecondWound = truncnorm::rtruncnorm(1, a = 0.5*30, mean = 2.5*30, sd = 2*30)) %>%
    ungroup() 
  
  # Outcomes for those in trajectory 8 (function of X1, X2, X4)
  trajectory8 <- temp %>%
    filter(trajectory == 8) %>%
    # Calculate time to death 
    rowwise() %>%
    mutate(timeToDeath = if_else(X4 == 0, 365*2*(log(runif(n = 1))/log(1-mort2yr_dom)),
                                 365*2*(log(runif(n = 1))/log(1-mort2yr_notDom)))) %>%
    # Calculate the time to amputation
    mutate(timeToAmp = if_else(X4 == 0, 365*2*(log(runif(n = 1))/log(1-amp2yr_dom)),
                               365*2*(log(runif(n = 1))/log(1-amp2yr_notDom)))) %>%
    # Calculate time to heal first wound
    mutate(timeToHealFirstWound = truncnorm::rtruncnorm(1, a = 0.5*30, mean = 2.5*30, sd = 2*30)) %>%
    # Calculate time to recurrence
    mutate(timeToRecur = if_else(X4 == 0, 365*2*(log(runif(n = 1))/log(1-recur2yr_dom)),
                                 365*2*(log(runif(n = 1))/log(1-recur2yr_notDom)))) %>%
    # Calculate time to heal second wound
    mutate(timeToHealSecondWound = truncnorm::rtruncnorm(1, a = 0.5*30, mean = 2.5*30, sd = 2*30)) %>%
    ungroup() 
  
  # Outcomes for those in trajectory 9 (function of X1)
  trajectory9 <- temp %>% 
    filter(trajectory == 9) %>%
    # Calculate time to death
    rowwise() %>%
    mutate(timeToDeath = if_else(X1 == 0, 365*2*(log(runif(n = 1))/log(1-mort2yr_dom)),
                                 365*2*(log(runif(n = 1))/log(1-mort2yr_notDom)))) %>%
    # Calculate time to amputation
    mutate(timeToAmp = if_else(X1 == 0, 365*2*(log(runif(n = 1))/log(1-amp2yr_dom)),
                               365*2*(log(runif(n = 1))/log(1-amp2yr_notDom)))) %>%
    # Calculate time to heal first wound
    mutate(timeToHealFirstWound = truncnorm::rtruncnorm(1, a = 0.5*30, mean = 2.5*30, sd = 2*30)) %>%
    ungroup() %>%
    mutate(timeToRecur = 10000) %>%
    mutate(timeToHealSecondWound = 0)
  
  # Outcomes for those in trajectory 10 (function of X1, X5)
  trajectory10 <- temp %>%
    filter(trajectory == 10) %>%
    # Calculate time to death 
    rowwise() %>%
    mutate(timeToDeath = if_else(X5 == 1, 365*2*(log(runif(n = 1))/log(1-mort2yr_dom)),
                                 365*2*(log(runif(n = 1))/log(1-mort2yr_notDom)))) %>%
    # Calculate the time to amputation
    mutate(timeToAmp = if_else(X5 == 1, 365*2*(log(runif(n = 1))/log(1-amp2yr_dom)),
                               365*2*(log(runif(n = 1))/log(1-amp2yr_notDom)))) %>%
    # Calculate time to heal first wound
    mutate(timeToHealFirstWound = truncnorm::rtruncnorm(1, a = 0.5*30, mean = 2.5*30, sd = 2*30)) %>%
    # Calculate time to recurrence
    mutate(timeToRecur = if_else(X5 == 1, 365*2*(log(runif(n = 1))/log(1-recur2yr_dom)),
                                 365*2*(log(runif(n = 1))/log(1-recur2yr_notDom)))) %>%
    # Calculate time to heal second wound
    mutate(timeToHealSecondWound = truncnorm::rtruncnorm(1, a = 0.5*30, mean = 2.5*30, sd = 2*30)) %>%
    ungroup() 
  
  # Outcomes for those in trajectory 11 (function of X1, X5)
  trajectory11 <- temp %>%
    filter(trajectory == 11) %>%
    # Calculate time to death 
    rowwise() %>%
    mutate(timeToDeath = if_else(X5 == 0, 365*2*(log(runif(n = 1))/log(1-mort2yr_dom)),
                                 365*2*(log(runif(n = 1))/log(1-mort2yr_notDom)))) %>%
    # Calculate the time to amputation
    mutate(timeToAmp = if_else(X5 == 0, 365*2*(log(runif(n = 1))/log(1-amp2yr_dom)),
                               365*2*(log(runif(n = 1))/log(1-amp2yr_notDom)))) %>%
    # Calculate time to heal first wound
    mutate(timeToHealFirstWound = truncnorm::rtruncnorm(1, a = 0.5*30, mean = 2.5*30, sd = 2*30)) %>%
    # Calculate time to recurrence
    mutate(timeToRecur = if_else(X5 == 0, 365*2*(log(runif(n = 1))/log(1-recur2yr_dom)),
                                 365*2*(log(runif(n = 1))/log(1-recur2yr_notDom)))) %>%
    # Calculate time to heal second wound
    mutate(timeToHealSecondWound = truncnorm::rtruncnorm(1, a = 0.5*30, mean = 2.5*30, sd = 2*30)) %>%
    ungroup() 
  

  
  # Combine get wound-free alive days and add censoring ------------------------
  trajectoriesWithTTE <- bind_rows(trajectory1, trajectory23, trajectory4, 
                                   trajectory5, trajectory6, trajectory7,
                                   trajectory8, trajectory9, trajectory10,
                                   trajectory11) %>% 
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