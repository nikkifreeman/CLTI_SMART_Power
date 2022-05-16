# Functions for Aim 1b power
# Author: Nikki Freeman
# Last modified: 4 May 2022

#' Calculate the expected pathway N
#'
#' @param studySize Integer; total number of participants in the study
#' @param designProbs vector of probabilities; Corresponds to the design probabilities
#' with designProbs = c(P_A, P_AA, P_AAA, P_ABA, P_BA)
#' @param responseProbs vector of probabilities; corresponds to the response probabilities
#' with responseProbs = c(P_1, P_2, P_3, P_4, P_5)--see the analysis plan for the definition 
#' of these probabilities; they also corresponds to each responder branch starting from top left 
#' and moving left to right and top to bottom
#'
#' @return vectore; vector of expected number of participants to go through 
#' each pathway given the simulation parameters; pathways are numbered from 1 to 11
#' working from the top to the bottom of the trial schematic
#' 
calculateExpectedPathwayN <- function(studySize, designProbs, responseProbs){
  n1 <- studySize * designProbs[1] * responseProbs[1]
  n2 <- studySize * designProbs[1] * (1-responseProbs[1]) * designProbs[2] * responseProbs[2]
  n3 <- studySize * designProbs[1] * (1-responseProbs[1]) * designProbs[2] * (1-responseProbs[2]) * responseProbs[3]
  n4 <- studySize * designProbs[1] * (1-responseProbs[1]) * designProbs[2] * (1-responseProbs[2]) * (1-responseProbs[3]) * designProbs[3]
  n5 <- studySize * designProbs[1] * (1-responseProbs[1]) * designProbs[2] * (1-responseProbs[2]) * (1-responseProbs[3]) * (1-designProbs[3])
  n6 <- studySize * designProbs[1] * (1-responseProbs[1]) * (1-designProbs[2]) * responseProbs[4]
  n7 <- studySize * designProbs[1] * (1-responseProbs[1]) * (1-designProbs[2]) * (1-responseProbs[4]) * designProbs[4]
  n8 <- studySize * designProbs[1] * (1-responseProbs[1]) * (1-designProbs[2]) * (1-responseProbs[4]) * (1-designProbs[4])
  n9 <- studySize * (1-designProbs[1]) * responseProbs[5]
  n10 <- studySize * (1-designProbs[1]) * (1-responseProbs[5]) * designProbs[5]
  n11 <- studySize * (1-designProbs[1]) * (1-responseProbs[5]) * (1-designProbs[5])
  out <- c(n1, n2, n3, n4, n5, n6, n7, n8, n9, n10, n11)
  names(out) <- paste0("path", 1:11)
  out
}


#' Calculate variance using the delta method
#'
#' @param delBeta vector; g'(x)
#' @param sandwichEst matrix; sandwich variance estimator
#'
#' @return numeric; variance 
getDeltaMethodVar <- function(delBeta, sandwichEst){
  delBeta %*% sandwichEst %*% delBeta
}

#' Calculate regime stadnard errors
#'
#' @param delBeta list; g'(x) for each  regime
#' @param mod model object; lm object from the regression
#'
#' @return vector; SEs for each regime
getRegimeSEs <- function(delBeta, mod){
  # Get the sandwhich estimator
  sandwichEst <- sandwich(mod)
  
  # Variances and standard errors
  regimeVariances <- map(delBeta, getDeltaMethodVar, sandwichEst = sandwichEst)
  regimeSEs <- sqrt(unlist(regimeVariances))
  
  regimeSEs
  
}


#' Add the weights to the observed data
#'
#' @param obsData dataframe; observed data as generated from generateSampleAim1b.R
#'
#' @return dataframe; contains obs data with `w` which is the likelihood of 
#' being assigned realized treatments and `w_inv` which is the inverse of w and 
#' the weight used in the regression
getObsDataWithWeights <- function(obsData){
  # Simplest censoring model b/c no covariates
  censoringModel <- glm(ltfu ~ 1, data = obsData)
  censoringWeights <- predict.glm(censoringModel, newData = obsData, type = "response")
  
  obsDataWithWeights <- obsData %>%
    mutate(censoringWeights = censoringWeights) %>%
    mutate(w = if_else(trajectory == 1, designProbs[1], -1),
           w = if_else(trajectory == 2, designProbs[1]*designProbs[2], w),
           w = if_else(trajectory == 3, designProbs[1]*designProbs[2], w),
           w = if_else(trajectory == 4, designProbs[1]*designProbs[2]*designProbs[3], w),
           w = if_else(trajectory == 5, designProbs[1]*designProbs[2]*(1-designProbs[3]), w),
           w = if_else(trajectory == 6, designProbs[1]*(1-designProbs[2]), w),
           w = if_else(trajectory == 7, designProbs[1]*(1-designProbs[2])*designProbs[4], w),
           w = if_else(trajectory == 8, designProbs[1]*(1-designProbs[2])*(1-designProbs[4]), w),
           w = if_else(trajectory == 9, (1-designProbs[1]), w),
           w = if_else(trajectory == 10, (1-designProbs[1])*designProbs[5], w),
           w = if_else(trajectory == 11, (1-designProbs[1])*(1-designProbs[5]), w)) %>%
    mutate(w = w * censoringWeights) %>%
    mutate(w_inv = 1/w)
  
  obsDataWithWeights
}


#' Create data with replications
#'
#' @param obsDataWithWeights dataframe; dataframe generated from getObsDataWithWeights()
#'
#' @return dataframe; contains the appropriate number of replications for each trajectory
getReplicatedData <- function(obsDataWithWeights){
  replicatedData <- bind_rows(filter(obsDataWithWeights, trajectory == 1) %>%
                                mutate(A2 = "A", A3 = "A"),
                              filter(obsDataWithWeights, trajectory == 1) %>%
                                mutate(A2 = "A", A3 = "D"),
                              filter(obsDataWithWeights, trajectory == 1) %>%
                                mutate(A2 = "B", A3 = "A"),
                              filter(obsDataWithWeights, trajectory == 1) %>%
                                mutate(A2 = "B", A3 = "D"),
                              filter(obsDataWithWeights, trajectory == 2) %>%
                                mutate(A3 = "A"),
                              filter(obsDataWithWeights, trajectory == 2) %>%
                                mutate(A3 == "A"),
                              filter(obsDataWithWeights, trajectory == 3) %>%
                                mutate(A3 = "A"),
                              filter(obsDataWithWeights, trajectory == 3) %>%
                                mutate(A3 = "D"),
                              filter(obsDataWithWeights, trajectory == 4),
                              filter(obsDataWithWeights, trajectory == 5),
                              filter(obsDataWithWeights, trajectory == 6) %>%
                                mutate(A3 == "A"),
                              filter(obsDataWithWeights, trajectory == 6) %>%
                                mutate(A3 == "D"),
                              filter(obsDataWithWeights, trajectory == 7),
                              filter(obsDataWithWeights, trajectory == 8),
                              filter(obsDataWithWeights, trajectory == 9) %>%
                                mutate(A2 == "A"),
                              filter(obsDataWithWeights, trajectory == 9) %>%
                                mutate(A3 == "D"),
                              filter(obsDataWithWeights, trajectory == 10),
                              filter(obsDataWithWeights, trajectory == 11))
  
  replicatedData
}

#' Get the analysis data
#'
#' @param replicatedData dataframe; as created from getReplicatedData()
#'
#' @return dataframe; dataframe with the indicators needed for the model
getAnalysisData <- function(replicatedData){
  analysisData <- replicatedData %>%
    # Create the Z vectors
    mutate(Z1 = if_else(A1 == "A", 1.0, 0.0),
           Z2 = if_else(A1 == "B", 1.0, 0.0),
           Z3 = if_else(A1 == "A" & A2 == "A", 1.0, 0.0),
           Z4 = if_else(A1 == "A" & A2 == "B", 1.0, 0.0)) %>%
    # Create the A10, A21, A22, A31, A32 vectors
    mutate(A10 = if_else(A1 == "A", 1, 0),
           A10 = if_else(A1 == "B", -1, A10)) %>%
    mutate(A21 = if_else(A1 == "A" & A2 == "A", 1, 0),
           A21 = if_else(A1 == "A" & A2 == "B", -1, A21)) %>%
    mutate(A22 = if_else(A1 == "B" & A2 == "A", 1, 0),
           A22 = if_else(A1 == "B" & A2 == "D", -1, A22)) %>%
    mutate(A31 = if_else(A1 == "A" & A2 == "A" & A3 == "A", 1, 0),
           A31 = if_else(A1 == "A" & A2 == "A" & A3 == "D", -1, A31)) %>%
    mutate(A32 = if_else(A1 == "A" & A2 == "B" & A3 == "A", 1, 0),
           A32 = if_else(A1 == "A" & A2 == "B" & A3 == "D", -1, A32))
  
  analysisData
}


#' Do one simulation run
#'
#' @param studySize numeric; total study size
#' @param designProbs vector of probabilities; from left to right and top down in 
#' the design schematic c(P_A, P_AA, P_AAA, P_ABA, P_BA)
#' @param responseProbRange list of vectors of range of responses; from left to right and top down
#' in the design schematic. Each vector has the form c(lower bound, upper bound) where the lower
#' bound is the lowest proportion of responders (w/ no recurrence) we'd expect and the upper bound
#' is the highest proportion of non-responders we'd expect
#' @param dominantRegime integer from 1 to 6; which regime should have the parameters
#' associated with "1" in the parameter grid; other regimes are associated with the 
#' "0" parameter values in the parameter grid
#' @param mort2yr_dom probability; 2 year mortality rate in the dominant regime
#' @param recur2yr_dom probability; 2 year recurrence rate in the dominant regime
#' @param mort2yr_notDom probability; 2 year mortality rate in the dominant regime
#' @param recur2yr_notDom probability; 2 year recurrence rate in the dominant regime
#' @param censorRate probability; ltfu rate
#'
#' @return
#' @export
#'
#' @examples
doOneSimRun <- function(studySize, designProbs, responseProbRange,
                        dominantRegime, mort2yr_dom,
                        recur2yr_dom, amp2yr_dom,
                        mort2yr_notDom, recur2yr_notDom, 
                        amp2yr_notDom, censorRate){
  # Generate simulated observed data -------------------------------------------
  obsData <- generateObserved(studySize = studySize, 
                              designProbs = designProbs, 
                              responseProbRange = responseProbRange,
                              dominantRegime = dominantRegime, 
                              mort2yr_dom = mort2yr_dom, 
                              recur2yr_dom = recur2yr_dom,
                              amp2yr_dom = amp2yr_dom,
                              mort2yr_notDom = mort2yr_notDom, 
                              recur2yr_notDom = recur2yr_notDom, 
                              amp2yr_notDom = amp2yr_notDom,
                              censorRate = censorRate)
  
  # Prepare for analysis -------------------------------------------------------
  # Calculate the relevant weights 
  obsDataWithWeights <- getObsDataWithWeights(obsData)
  
  # Replicate 
  replicatedData <- getReplicatedData(obsDataWithWeights)
  
  # Create analysis data (make indicators)
  analysisData <- getAnalysisData(replicatedData = replicatedData) 
  
  # Estimate model -------------------------------------------------------------
  mod <- lm(woundFreeAliveDays ~ A10 + Z1:A21 +Z2:A22 + Z3:A31 + Z4:A32, 
            weights = w_inv, data = analysisData)
  
  # Analyze estimated model ----------------------------------------------------
  # Regime means
  ## Regime coefficients
  delBeta <- list(delBeta1 = c(1, 1, 1, 0, 1, 0), 
                  delBeta2 = c(1, 1, 1, 0, -1, 0),
                  delBeta3 = c(1, 1, -1, 0, 0, 1),
                  delBeta4 = c(1, 1, -1, 0, 0, -1),
                  delBeta5 = c(1, -1, 0, 1, 0, 0),
                  delBeta6 = c(1, -1, 0, -1, 0, 0))
  # Estimated value of embedded regimes
  regimeMeans <- data.frame(matrix(unlist(map(delBeta, ~ sum(. * summary(mod)$coefficients[,1]))), nrow = 1))
  names(regimeMeans) <- paste0("regimeMean", 1:6)
  
  # P-value (using F-test, somewhat robust to non-normality if N large)
  fstat <- summary(mod)$fstatistic
  pval <- data.frame(pvalue = pf(fstat[1], fstat[2], fstat[3], lower.tail = FALSE))
  
  # Estimated SE of embedded regimes
  regimeSEs <- data.frame(matrix(getRegimeSEs(delBeta = delBeta, mod = mod), nrow = 1))
  names(regimeSEs) <- paste0("regimeSE", 1:6)
  
  # Things to return: regimeMeans, regimeSEs, pval
  n_trajectory <- obsData %>% group_by(trajectory) %>% count() %>%
    mutate(trajectory = paste0('n_trajectory', trajectory)) %>%
    pivot_wider(names_from = trajectory, values_from = n)
  
  bind_cols(regimeMeans, regimeSEs, pval, n_trajectory)
}

