# Generate a simulation sample --including parameterization for amputation rate
# Author: Nikki Freeman
# Last modified: 9 May 2022

# Preamble ---------------------------------------------------------------------
# These are the functions for generating a single sample for the simulation
# analysis. 
# Required packages: tidyverse
# Required scripts: 1_generateSampleAim1a.R

# Function ---------------------------------------------------------------------

#' Generate a single simulation data set for Aim 1.1 that includes amputation
#'
#' @param N positive integer; same N used for baseline data, sample size
#' @param p_trt1 real in (0, 1); probability of recieving "treatment 1"
#' @param mort0 real in (0, 1); Pr(2 year mortality) in arm 0
#' @param mort1 real in (0, 1); Pr(2 year mortality) in arm 1
#' @param recur0 real in (0, 1); Pr(recurrence in 2 years) in arm 0
#' @param recur1 real in (0, 1); Pr(recurrence in 2 years) in arm 1
#' @param amp0 real in (0, 1); Pr(amputation in 2 years) in arm 0
#' @param amp1 real in (0, 1); Pr(amputation in 2 years) in arm 1
#' @param C_L real in (0, 1); Pr(ltfu over the study period)
#' @param beta_mort real[4]; coefficients for X1:X4 in the TTE mortality 
#' @param beta_recur real[4]; coefficients for X1:X4 in the TTE recurrence
#' @param beta_amp real[4]; coefficients for X1:X4 in the TTE amputation
#'
#' @return dataframe; contains 1 simulation data set for aim 1.1
simulateDataAim1aWithAmp_exponential <- function(N, p_trt1, 
                                          mort0, mort1, 
                                          recur0, recur1, 
                                          amp0, amp1, C_L, 
                                          beta_mort, beta_recur, beta_amp){
  ltfu_ids <- sample(1:N, size = round(N*C_L), replace = FALSE)
  temp <- generateBaselineData(N) %>%
    assignFirstRandomization(., N = N, p_trt1 = p_trt1) %>%
    mutate(linPred_mort = beta_mort[1]*X1 + beta_mort[2]*X2 + beta_mort[3]*X3 + beta_mort[4]*X4,
           expLinPred_mort = exp(linPred_mort)) %>%
    mutate(linPred_recur = beta_recur[1]*X1 + beta_recur[2]*X2 + beta_recur[3]*X3 + beta_recur[4]*X4,
           expLinPred_recur = exp(linPred_recur)) %>%
    mutate(linPred_amp = beta_amp[1]*X1 + beta_amp[2]*X2 + beta_amp[3]*X3 + beta_amp[4]*X4,
           expLinPred_amp = exp(linPred_amp)) %>%
    rowwise() %>%
    mutate(timeToDeath =
             if_else(A1 == 0,
                     generateEventTime_exponential(mort0, t = 365*2, expLinPred_mort),
                     generateEventTime_exponential(mort1, t = 365*2, expLinPred_mort))) %>%
    mutate(timeToAmp = 
             if_else(A1 == 0, 
                     generateEventTime_exponential(0.2, t = 365*2, expLinPred_amp), 
                     generateEventTime_exponential(0.25, t = 365*2, expLinPred_amp))) %>%
    mutate(timeToRecur =
             if_else(A1 == 0,
                     generateEventTime_exponential(recur0, t = 365*2, expLinPred_recur),
                     generateEventTime_exponential(recur1, t = 365*2, expLinPred_recur))) %>%
    mutate(timeToHealFirstWound = truncnorm::rtruncnorm(1, a = 0.5*30, mean = 4.5*30, sd = 2*30)) %>%
    mutate(timeToHealSecondWound = truncnorm::rtruncnorm(1, a = 0.5*30, mean = 4.5*30, sd = 2*30)) %>%
    ungroup() %>%
    mutate(ltfu = if_else(id %in% ltfu_ids, 1, 0)) %>%
    mutate(firstWoundHealingDay = timeToHealFirstWound) %>%
    rowwise() %>%
    mutate(timeToRecur = max(timeToRecur, firstWoundHealingDay + 1)) %>%
    ungroup() %>%
    mutate(secondWoundHealingDay = timeToRecur + timeToHealSecondWound) %>%
    mutate(endOfStudy = 365*2)

  # for convenience, make event1:event5 variables
  sequences <- temp %>%
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


  allDf <- left_join(temp, sequences, by = "id")

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
