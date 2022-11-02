# 2_aim2PowerFnxs.R
# Author: Nikki Freeman
# Last modified: 20 June 2022

# Packages ---------------------------------------------------------------------
library(tidyverse)

# Scripts ----------------------------------------------------------------------
if(str_detect(getwd(), "resubmission")){
  source("./3_Aim2/1_code/1_generateSimulationData.R")
} else{
  source("1_generateSimulationData.R")
}

out <- data.frame(k = vector(), value = numeric())

tic()
for(l in 1:1000){
  obsData <- generateObserved(studySize = 650, 
                              designProbs = rep(0.5, 5), 
                              responseProbRange = list(p1 = c(0.25, 0.35),
                                                       p2 = c(0.25, 0.35),
                                                       p3 = c(0.3, 0.6),
                                                       p4 = c(0.3, 0.6),
                                                       p5 = c(0.3, 0.6)), 
                              dominantRegime = 5, 
                              mort2yr_dom = 0.15, recur2yr_dom = 0.3, amp2yr_dom = 0.2,
                              mort2yr_notDom = 0.20, recur2yr_notDom = 0.35, amp2yr_notDom = 0.25, 
                              censorRate = 0.15)
  
  # Create the training data and the evaluation data -----------------------------
  # Create the folds
  nFolds <- 2
  folds_list <- caret::createFolds(obsData$woundFreeAliveDays, k = nFolds)
  
  for(i in 1:nFolds){
    trainingData <- obsData %>% filter(!(id %in% folds_list[[i]]))
    testingData <- obsData %>% filter(id %in% folds_list[[i]])
    
    # Learn the optimal policy ---------------------------------------------------
    # Stage 3 Q-learning for upper randomization -----------------------------------
    # Data for the 3rd stage upper randomization
    dat_stage3upper <- trainingData %>% filter(trajectory == 4 | trajectory == 5) %>%
      mutate(A3 = if_else(A3 == "A", 1, 0)) %>%
      mutate(X1A3 = X1*A3, 
             X2A3 = X2*A3,
             X3A3 = X3*A3,
             X4A3 = X4*A3,
             X5A3 = X5*A3)
    # Estimate the Q-function
    mod_stage3upper <- lm(woundFreeAliveDays ~ X1 + X2 + X3 + X4 + X5 +  A3 + 
                            X1A3 + X2A3 + X3A3 + X4A3 + X5A3, 
                          data = dat_stage3upper)
    # Prediction for A3 == 1 (A3 == A)
    dat_stage3upper <- trainingData %>% filter(trajectory == 4 | trajectory == 5) %>%
      mutate(A3 = 1) %>%
      mutate(X1A3 = X1*A3, 
             X2A3 = X2*A3,
             X3A3 = X3*A3,
             X4A3 = X4*A3,
             X5A3 = X5*A3)
    pred_stage3upperA <- predict.lm(mod_stage3upper, dat_stage3upper)
    # Prediction for A3 == 0 (A3 == D)
    dat_stage3upper <- trainingData %>% filter(trajectory == 4 | trajectory == 5) %>%
      mutate(A3 = 0) %>%
      mutate(X1A3 = X1*A3, 
             X2A3 = X2*A3,
             X3A3 = X3*A3,
             X4A3 = X4*A3,
             X5A3 = X5*A3)
    pred_stage3upperD <- predict.lm(mod_stage3upper, newdata = dat_stage3upper)
    
    # Replace wound free alive days with pseudo values
    dat_stage3upper$woundFreeAliveDays[pred_stage3upperA >= pred_stage3upperD] <- pred_stage3upperA[pred_stage3upperA >= pred_stage3upperD]
    dat_stage3upper$woundFreeAliveDays[pred_stage3upperA < pred_stage3upperD] <- pred_stage3upperD[pred_stage3upperA < pred_stage3upperD]
    
    # data set with pseudo values
    updatedData <- left_join(trainingData, dat_stage3upper %>% select(id, pseudo = woundFreeAliveDays), 
                             by = "id") %>% 
      mutate(woundFreeAliveDays = if_else(!is.na(pseudo), pseudo, woundFreeAliveDays)) %>%
      select(-pseudo)
    
    # Stage 3 Q_learning for lower randomization -----------------------------------
    # Data for the 3rd stage lower randomization
    dat_stage3lower <- trainingData %>% filter(trajectory == 7 | trajectory == 8) %>%
      mutate(A3 = if_else(A3 == "A", 1, 0)) %>%
      mutate(X1A3 = X1*A3, 
             X2A3 = X2*A3,
             X3A3 = X3*A3,
             X4A3 = X4*A3,
             X5A3 = X5*A3)
    # Estimate the Q-function
    mod_stage3lower <- lm(woundFreeAliveDays ~ X1 + X2 + X3 + X4 + X5 +  A3 + 
                            X1A3 + X2A3 + X3A3 + X4A3 + X5A3, 
                          data = dat_stage3lower)
    # Prediction for A3 == 1 (A3 == A)
    dat_stage3lower <- trainingData %>% filter(trajectory == 7 | trajectory == 8) %>%
      mutate(A3 = 1) %>%
      mutate(X1A3 = X1*A3, 
             X2A3 = X2*A3,
             X3A3 = X3*A3,
             X4A3 = X4*A3,
             X5A3 = X5*A3)
    pred_stage3lowerA <- predict.lm(mod_stage3lower, dat_stage3lower)
    # Prediction for A3 == 0 (A3 == D)
    dat_stage3lower <- trainingData %>% filter(trajectory == 7 | trajectory == 8) %>%
      mutate(A3 = 0) %>%
      mutate(X1A3 = X1*A3, 
             X2A3 = X2*A3,
             X3A3 = X3*A3,
             X4A3 = X4*A3,
             X5A3 = X5*A3)
    pred_stage3lowerD <- predict.lm(mod_stage3lower, newdata = dat_stage3lower)
    
    # Replace wound free alive days with pseudo values
    dat_stage3lower$woundFreeAliveDays[pred_stage3lowerA >= pred_stage3lowerD] <- pred_stage3lowerA[pred_stage3lowerA >= pred_stage3lowerD]
    dat_stage3lower$woundFreeAliveDays[pred_stage3lowerA < pred_stage3lowerD] <- pred_stage3lowerD[pred_stage3lowerA < pred_stage3lowerD]
    
    # data set with pseudo values
    updatedData <- left_join(updatedData, dat_stage3lower %>% select(id, pseudo = woundFreeAliveDays), 
                             by = "id") %>% 
      mutate(woundFreeAliveDays = if_else(!is.na(pseudo), pseudo, woundFreeAliveDays)) %>%
      select(-pseudo)
    
    # Stage 2 Q-learning for the upper randomization -------------------------------
    # Data for the 2nd stage upper randomization
    dat_stage2upper <- updatedData %>% filter(trajectory %in% c(2, 3, 4, 5, 6, 7, 8)) %>%
      mutate(A2 = if_else(A2 == "A", 1, 0)) %>%
      mutate(X1A2 = X1*A2, 
             X2A2 = X2*A2,
             X3A2 = X3*A2,
             X4A2 = X4*A2,
             X5A2 = X5*A2)
    # Estimate the Q-function
    mod_stage2upper <- lm(woundFreeAliveDays ~ X1 + X2 + X3 + X4 + X5 +  A2 + 
                            X1A2 + X2A2 + X3A2 + X4A2 + X5A2, 
                          data = dat_stage2upper)
    # Prediction for A2 == 1 (A2 == A)
    dat_stage2upper <- updatedData %>% filter(trajectory %in% c(2, 3, 4, 5, 6, 7, 8)) %>%
      mutate(A2 = 1) %>%
      mutate(X1A2 = X1*A2, 
             X2A2 = X2*A2,
             X3A2 = X3*A2,
             X4A2 = X4*A2,
             X5A2 = X5*A2)
    pred_stage2upperA <- predict.lm(mod_stage2upper, dat_stage2upper)
    # Prediction for A3 == 0 (A3 == B)
    dat_stage2upper <- updatedData %>% filter(trajectory %in% c(2, 3, 4, 5, 6, 7, 8)) %>%
      mutate(A2 = 0) %>%
      mutate(X1A2 = X1*A2, 
             X2A2 = X2*A2,
             X3A2 = X3*A2,
             X4A2 = X4*A2,
             X5A2 = X5*A2)
    pred_stage2upperB <- predict.lm(mod_stage2upper, newdata = dat_stage2upper)
    
    # Replace wound free alive days with pseudo values
    dat_stage2upper$woundFreeAliveDays[pred_stage2upperA >= pred_stage2upperB] <- pred_stage2upperA[pred_stage2upperA >= pred_stage2upperB]
    dat_stage2upper$woundFreeAliveDays[pred_stage2upperA < pred_stage2upperB] <- pred_stage2upperB[pred_stage2upperA < pred_stage2upperB]
    
    # data set with pseudo values
    updatedData <- left_join(updatedData, dat_stage2upper %>% select(id, pseudo = woundFreeAliveDays), 
                             by = "id") %>% 
      mutate(woundFreeAliveDays = if_else(!is.na(pseudo), pseudo, woundFreeAliveDays)) %>%
      select(-pseudo)
    
    # Stage 2 Q-learning for the lower randomization -------------------------------
    # Data for the 2nd stage lower randomization
    dat_stage2lower <- updatedData %>% filter(trajectory %in% c(10, 11)) %>%
      mutate(A2 = if_else(A2 == "A", 1, 0)) %>%
      mutate(X1A2 = X1*A2, 
             X2A2 = X2*A2,
             X3A2 = X3*A2,
             X4A2 = X4*A2,
             X5A2 = X5*A2)
    # Estimate the Q-function
    mod_stage2lower <- lm(woundFreeAliveDays ~ X1 + X2 + X3 + X4 + X5 +  A2 + 
                            X1A2 + X2A2 + X3A2 + X4A2 + X5A2, 
                          data = dat_stage2lower)
    # Prediction for A2 == 1 (A2 == A)
    dat_stage2lower <- updatedData %>% filter(trajectory %in% c(10, 11)) %>%
      mutate(A2 = 1) %>%
      mutate(X1A2 = X1*A2, 
             X2A2 = X2*A2,
             X3A2 = X3*A2,
             X4A2 = X4*A2,
             X5A2 = X5*A2)
    pred_stage2lowerA <- predict.lm(mod_stage2lower, dat_stage2lower)
    # Prediction for A3 == 0 (A3 == D)
    dat_stage2lower <- updatedData %>% filter(trajectory %in% c(10, 11)) %>%
      mutate(A2 = 0) %>%
      mutate(X1A2 = X1*A2, 
             X2A2 = X2*A2,
             X3A2 = X3*A2,
             X4A2 = X4*A2,
             X5A2 = X5*A2)
    pred_stage2lowerD <- predict.lm(mod_stage2lower, newdata = dat_stage2lower)
    
    # Replace wound free alive days with pseudo values
    dat_stage2lower$woundFreeAliveDays[pred_stage2lowerA >= pred_stage2lowerD] <- pred_stage2lowerA[pred_stage2lowerA >= pred_stage2lowerD]
    dat_stage2lower$woundFreeAliveDays[pred_stage2lowerA < pred_stage2lowerD] <- pred_stage2lowerD[pred_stage2lowerA < pred_stage2lowerD]
    
    # data set with pseudo values
    updatedData <- left_join(updatedData, dat_stage2lower %>% select(id, pseudo = woundFreeAliveDays), 
                             by = "id") %>% 
      mutate(woundFreeAliveDays = if_else(!is.na(pseudo), pseudo, woundFreeAliveDays)) %>%
      select(-pseudo)
    
    # Stage 1 Q-learning -----------------------------------------------------------
    # Data for the 1st stage randomization
    dat_stage1 <- updatedData %>% 
      mutate(A1 = if_else(A1 == "A", 1, 0)) %>%
      mutate(X1A1 = X1*A1, 
             X2A1 = X2*A1,
             X3A1 = X3*A1,
             X4A1 = X4*A1,
             X5A1 = X5*A1)
    # Estimate the Q-function
    mod_stage1 <- lm(woundFreeAliveDays ~ X1 + X2 + X3 + X4 + X5 +  A1 + 
                       X1A1 + X2A1 + X3A1 + X4A1 + X5A1, 
                     data = dat_stage1)
    # Prediction for A2 == 1 (A2 == A)
    dat_stage1 <- updatedData %>% 
      mutate(A1 = 1) %>%
      mutate(X1A1 = X1*A1, 
             X2A1 = X2*A1,
             X3A1 = X3*A1,
             X4A1 = X4*A1,
             X5A1 = X5*A1)
    pred_stage1A <- predict.lm(mod_stage1, dat_stage1)
    # Prediction for A3 == 0 (A3 == D)
    dat_stage1 <- updatedData %>%
      mutate(A1 = 0) %>%
      mutate(X1A1 = X1*A1, 
             X2A1 = X2*A1,
             X3A1 = X3*A1,
             X4A1 = X4*A1,
             X5A1 = X5*A1)
    pred_stage1B <- predict.lm(mod_stage1, newdata = dat_stage1)
    
    # Replace wound free alive days with pseudo values
    dat_stage1$woundFreeAliveDays[pred_stage1A >= pred_stage1B] <- pred_stage1A[pred_stage1A >= pred_stage1B]
    dat_stage1$woundFreeAliveDays[pred_stage1A < pred_stage1B] <- pred_stage1B[pred_stage1A < pred_stage1B]
    
    # data set with pseudo values
    updatedData <- left_join(updatedData, dat_stage1 %>% select(id, pseudo = woundFreeAliveDays), 
                             by = "id") %>% 
      mutate(woundFreeAliveDays = if_else(!is.na(pseudo), pseudo, woundFreeAliveDays)) %>%
      select(-pseudo)
    
    # Evaluate the learned policy on the test set
    baselineCovariateMatrix <- testingData %>% 
      mutate(intercept = 1) %>% 
      mutate(A1 = if_else(A1 == "A", 1, 0),
             A2 = if_else(A2 == "A", 1, 0),
             A3 = if_else(A3 == "A", 1, 0)) %>%
      select(intercept, X1:X5, A1, A2, A3)
    
    # Optimal rule at upper decision 3
    covariateMatrix3 <-baselineCovariateMatrix %>% select(-c(A1, A2)) %>%
      mutate(X1A3 = X1*A3,
             X2A3 = X2*A3,
             X3A3 = X3*A3,
             X4A3 = X4*A3,
             X5A3 = X5*A3) %>%
      select(-c(intercept, X1, X2, X3, X4, X5))
    covariateMatrix3 <- as.matrix(covariateMatrix3)
    A3upper_opt <- (covariateMatrix3 %*% mod_stage3upper$coefficients[7:12] > 0)*1
    
    # Optimal rule at lower decision 3
    A3lower_opt <- (covariateMatrix3 %*% mod_stage3lower$coefficients[7:12] > 0)*1
    
    # Optimal rule at upper decision 2
    covariateMatrix2 <- baselineCovariateMatrix %>%
      select(-c(A1, A3)) %>%
      mutate(X1A2 = X1*A2,
             X2A2 = X2*A2,
             X3A2 = X3*A2,
             X4A2 = X4*A2,
             X5A2 = X5*A2) %>%
      select(-c(intercept, X1, X2, X3, X4, X5))
    covariateMatrix2 <- as.matrix(covariateMatrix2)
    A2upper_opt <- (covariateMatrix2 %*% mod_stage2upper$coefficients[7:12] >0 )*1
    
    # Optimal rule at lower decision 2
    A2lower_opt <- (covariateMatrix2 %*% mod_stage2upper$coefficients[7:12] > 0)*1
    
    # Optimal rule at decision 1
    covariateMatrix1 <- baselineCovariateMatrix %>%
      select(-c(A2, A3)) %>%
      mutate(X1A1 = X1*A1,
             X2A1 = X2*A1,
             X3A1 = X3*A1,
             X4A1 = X4*A1,
             X5A1 = X5*A1) %>%
      select(-c(intercept, X1, X2, X3, X4, X5))
    covariateMatrix1 <- as.matrix(covariateMatrix1)
    A1_opt <- (covariateMatrix1 %*% mod_stage1$coefficients[7:12] > 0)*1
    
    valueDf <- testingData %>%
      mutate(A1Opt = A1_opt,
             A2lowerOpt = A2lower_opt,
             A2upperOpt = A2upper_opt,
             A3lowerOpt = A3lower_opt,
             A3upperOpt = A3upper_opt) %>%
      mutate(consistent1 = if_else(A1Opt == 1 & A1 == "A", 1, 0)) %>%
      mutate(consistent2 = if_else(A1Opt == 0 & A1 == "B", 1, consistent1)) %>%
      mutate(consistent2 = if_else(trajectory == 1, 1, 0)) %>%
      mutate(consistent3 = if_else(trajectory == 1, 1, 0)) %>%
      mutate(consistent2 = if_else((trajectory == 23) & A2upperOpt == 1 & A2 == "A", 1, consistent2)) %>%
      mutate(consistent3 = if_else(trajectory == 23, 1, consistent3)) %>%
      mutate(consistent2 = if_else(trajectory == 4 & A2upperOpt == 1 & A2 == "A", 1, consistent2)) %>%
      mutate(consistent3 = if_else(trajectory == 4 & A3upperOpt == 1 & A3 == "A", 1, consistent3)) %>%
      mutate(consistent2 = if_else(trajectory == 5 & A2upperOpt == 1 & A2 == "A", 1, consistent2)) %>%
      mutate(consistent3 = if_else(trajectory == 5 & A3upperOpt == 0 & A3 == "D", 1, consistent3)) %>%
      mutate(consistent2 = if_else(trajectory == 6 & A2upperOpt == 0 & A2 == "B", 1, consistent2)) %>%
      mutate(consistent3 = if_else(trajectory == 6, 1, consistent3)) %>%
      mutate(consistent2 = if_else(trajectory == 7 & A2upperOpt == 0 & A2 == "B", 1, consistent2)) %>%
      mutate(consistent3 = if_else(trajectory == 7 & A3lowerOpt == 1 & A3 == "A", 1, consistent3)) %>%
      mutate(consistent2 = if_else(trajectory == 8 & A2upperOpt == 0 & A2 == "B", 1, consistent2)) %>%
      mutate(consistent3 = if_else(trajectory == 8 & A3lowerOpt == 0 & A3 == "D", 1, consistent3)) %>%
      mutate(consistent2 = if_else(trajectory == 9, 1, consistent2)) %>%
      mutate(consistent3 = if_else(trajectory == 9, 1, consistent3)) %>%
      mutate(consistent2 = if_else(trajectory == 10 & A2lowerOpt == 1 & A2 == "A", 1, consistent2)) %>%
      mutate(consistent3 = if_else(trajectory == 10, 1, consistent3)) %>%
      mutate(consistent2 = if_else(trajectory == 11 & A2lowerOpt == 0 & A2 == "D", 1, consistent2)) %>%
      mutate(consistent3 = if_else(trajectory == 11, 1, consistent3)) %>%
      mutate(weight = if_else(trajectory %in% c(1, 9), 2, 0),
             weight = if_else(trajectory %in% c(23, 6, 10, 11), 4, weight),
             weight = if_else(trajectory %in% c(4, 5, 7, 8), 8, weight)) %>%
      mutate(summand = weight*woundFreeAliveDays*consistent1*consistent2*consistent3) %>%
      mutate(stabilization = (consistent1*consistent2*consistent3)*weight) %>%
      summarise(value = sum(summand)/sum(stabilization)) # Value = 552.4672, 547.1429 (avg = 549.8051)
    
    out <- bind_rows(out, data.frame(l = l, k = i, value = valueDf$value, mort2yr_dom = 0.15, recur2yr_dom = 0.3, amp2yr_dom = 0.2,
                                     mort2yr_notDom = 0.20, recur2yr_notDom = 0.35, amp2yr_notDom = 0.25, 
                                     censorRate = 0.15))
  }
  
}

toc()

write_csv(out, file = "./3_Aim2/2_pipeline/2_aim2Power_1.csv")


