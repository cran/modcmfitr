
# Lines for development only
#   rm(list=ls(all=TRUE))
#   library(nloptr)
#   library(gtools)
# End development lines

# Quantile of Scaled Beta
#
# This function returns the quantile of a scaled beta distribution.  A scaled beta is a standard beta distribution, rescaled between the desired lower and upper limits.
# @param x Vector of 5 parameters: c(quantile, alpha, beta, lower limit, upper limit)
# @keywords scaled beta
# @examples
# qScaledBeta(c(0.025,1,2,0.5,1))
qScaledBeta <- function(x) {
  q <- x[4]+qbeta(x[1],x[2],x[3])*(x[5]-x[4])
  return(q)
}

# Sample from Scaled Beta distribution
#
# This function returns n samples from a scaled beta. A scaled beta is a standard beta distribution, rescaled between the desired lower and upper limits.
# @param x Vector of 5 parameters: c(n, alpha, beta, lower limit, upper limit)
# @keywords scaled beta
# @examples
# rScaledBeta(c(1000,1,2,0.5,1))
rScaledBeta <- function(x) {
  x<- as.numeric(x)
  r <- rbeta(x[1],x[2],x[3])
  r <- x[4]+r*(x[5]-x[4])
  return(r)
}

#' Sample from modified Connor-Mosimann distribution
#'
#' This function returns n samples from a modified CM distribution
#' @param n Number of samples
#' @param Z vector of (k-1) scaled beta parameters.  Where k=number of dimensions to problem. Columns = (a,b,lower limit, upper limit) of scaled beta distribution
#' @keywords modified Connor-Mosimann
#' @export
#' @examples
#' n <- 1000
#' Z <- matrix(data = c(1.098, 3.070, 0.66, 0.85,
#'                    9.021, 7.990, 0.62, 1,
#'                    10, 0.01, 0.14, 0.58,
#'                    10, 10, 0.57, 1),ncol=4,byrow=TRUE)
#' rModCM(n,Z)
rModCM <- function(n,Z) {
  #solve bug whereby if Z only has one row then it is imported as a list not a matrix
  #Z <- matrix(Z,ncol=4)

  # count how many dimensions:
  ewDims <- nrow(Z)+1

  #solve error whereby single row does not return valid nrow()
  if (length(ewDims)==0) ewDims <- 2
##################################################################################

  #create matrices to hold sampled Zs and Ps
  Zval <- matrix(data=NA,nrow=n,ncol=ewDims)
  Pval <- matrix(data=NA,nrow=n,ncol=ewDims)

  if (ewDims == 2) {
    Zval[,1] <- rScaledBeta(c(n,Z))
  } else {
    for (i in 1:(ewDims-1)) {
      Zval[,i] <- rScaledBeta(c(n,Z[i,]))
    }
  }
  Zval[,ewDims] <- 1

  Pval[,1] <- Zval[,1]
  for (i in 2:ewDims) {
    scaling <- 1
    for (j in (i-1):1) {
      scaling <- scaling*(1-Zval[,j])
    }
    Pval[,i] <- Zval[,i]*scaling
  }
  return(Pval)
}

# Check beta distributions are consistent
#
# This function inverts parameters of a scaled beta distribution if the lower limit is above the upper limit.
# This is needed because the search algorithm used in fitModCM() does not impose restrictions on LL<=UL.
# @param x Vector of scaled beta hyperparameters: (a,b,lower limit, upper limit)
# @keywords scaled beta
# @examples
#
# consistentBetas(c(1,2,0.8,0.6))
consistentBetas <-function(x) {
  ifelse (x[3]>x[4], return(c(x[2],x[1],x[4],x[3])), return(x))
}

# Calculate sum of squared deviation between modelled and elicited quantiles for modified Connor-Mosimann distribution
#
# This function uses a set of hyperparameters of an mCM distribution as an input and calculates the
# associated quantiles of the marginal scaled betas for that distribution, then compares this
# with the target quantiles by calculating the sum of squared deviation.
# The function returns the SSD.
# @param x Vector of input parameters: NrDims, Iterations, Quantiles being compared,
# target for quantiles and hyperparameters of mCM distribution.
# Note input must be a single atomic vector not dataframes or matrices.  This is because the
# search algorithm crs2lm(), called in fitmodCM() can only send a single vector when it calls this function.
# Parameters of Zed to try out (a1,b1,L1,U1)... etc
# @param NrDims Number of dimensions to mCM distribution (e.g. number of outcomes such as remission, progression, dead)
# @param Iterations Number of Samples from the mCM distribution to pull to calculate the quantiles.  The higher the more precise the reported quantiles.
# @param Quantiles Vector of quantiles to be matching.
# @param Targets Target values for lower, middle and upper quantiles for each dimension (i.e. quantities elicited from experts)
# @param Zeds Hyperparamters of mCM distribution from which to calculate the quantiles.
# @keywords modified Connor-Mosimann
# @examples
# NrDims <- 3
# Iterations <-10000
# Quantiles <- c(0.025,0.5,0.975) #(median and upper and lower 95% credibility interval)
# Targets <- c(0.3,0.4,0.5,
#              0.1,0.4,0.5,
#              0.05,0.2,0.3)
# Zeds <- c(1,1,0,1,
#          1,1,0,1,
#          1,1,0,1)
# TestFitModCM(c(NrDims,Iterations,Quantiles,Targets,Zeds))
#
# 
TestFitModCM <- function(x) {
  x <- as.numeric(x)
  #assemble input vector into something easier to understand
  NrDims <- x[1]
  Iterations <- x[2]
  Quantiles <- x[3:5]
  Targets <- matrix(data=x[6:(5+3*NrDims)], nrow = NrDims, ncol = 3, byrow = TRUE)
  Zeds <- matrix(data=x[(6+3*NrDims):(length(x))], nrow = NrDims, ncol = 4, byrow = TRUE)

  #sample from Modified CM dbn
  p<-rModCM(Iterations,Zeds)

  #Find the quantiles from the sampled values
  ewResults <- matrix(quantile(p[,1], Quantiles), nrow = NrDims, ncol = 3, byrow=TRUE)
  for (i in 2:NrDims) {
    ewResults[i,] <- quantile(p[,i], Quantiles)
  }

  #Code to place greater weight on getting the medians right over the LL & ULs.  To weight uncomment first line.
  #to leave unweighted, uncomment second.
  #ewWeights <- matrix(rep(c(1,NrDims,1),NrDims), nrow = NrDims, ncol = 3, byrow=TRUE)
  ewWeights <- matrix(rep(c(1,1,1),NrDims), nrow = NrDims, ncol = 3, byrow=TRUE)


  #calculate (weighted) sum of the squared deviation from target
  SSD = sum((ewWeights*(Targets-ewResults))^2)
  cat("\n",SSD)
  return(SSD)
}

#' Fit modified Connor-Mosimann or Connor-Mosimann distribution
#'
#' This function fits a modified CM distribution to an elicited set of quantiles from, e.g. an expert elicitation workshop.
#' The function uses the crs2lm() function from the nloptr package to search for the
#' set of hyperparameters that that generates quantiles that match the elicited data as closely as possible.
#' crs2lm() repeatedly calls TestFitModCM() until it finds the best fitting set of inputs.
#' The function returns a data frame containing the best fit set of hyperparameters for the Zeds, the Target quantiles and the modelled quantiles.
#'
#' @param Outcomes A vector of names of outcomes
#' @param RawData matrix of lower, middle and upper quantiles for each dimension elicited from experts
#' @param SearchParams A vector of number of iterations and max number of searches.
#' Number of iterations is the number of draws from the mCM distribution used to estimate the quantiles.  Try 10,000 first.
#' Max number of searches is the maximum number of searches the search algorithm conducts.  The higher the better the fit, but also the longer it takes to find it.  
#' @param ModCMorCM If =1, then will fit a modified CM.  If =0 then will fit a CM (i.e. forces lower and upper limits of each marginal scaled beta to 0 and 1 respectively)
#' @param Quantiles Sets the quantiles to be fit.  If median and 95\% Credibility Interval, then set to c(0.025,0.5,0.975).  If median and tertiles then (0.33,0.5,0.66).  If median and quartiles then c(0.25,0.5,0.75) and so on.
#' @keywords modified Connor-Mosimann
#' @export
#' @examples
#' Outcomes <- c("Remission","Progression","Dead")
#' RawData <- matrix(data = c(0.43, 0.55, 0.65,
#'                           0.16, 0.27, 0.46,
#'                           0.03, 0.18, 0.23
#'            ),ncol=3,byrow=TRUE)
#'
#' SearchParams <- c(10000,100) #number of iterations, max number of searches
#'
#' ModCMorCM <- 1
#'
#' Quantiles <- c(0.025,0.5,0.975) # example here is 95% credibility limits and median.
#'
#' fitModCM(Outcomes, RawData, SearchParams, ModCMorCM, Quantiles)
#' @return Returns a matrix, each row representing one of the Outcomes.  Columns a, b, LL and UL are the parameters of the mCM distribution.  
#' Note that the final row will always be zeroes as the number of rows of the mCM distribution is k-1 where k is the number of outcomes.
#' The next three columns (Tgt_LL, Tgt_MED, Tgt_UL) are the target quantiles input to the function as 'RawData' in the example above.
#' The final three (Mdl_LL, Mdl_MED, Mdl_UL) are the quantiles resulting from the model fitting.  If the mCM is a good fit, the Mdl columns will be 
#' identical to the Tgt columns. 
fitModCM <- function(Outcomes, RawData, SearchParams, ModCMorCM, Quantiles) {

  NrDims <- length(Outcomes)
  Iterations <- SearchParams[1]
  MaxSearches <- SearchParams[2]

  #If statement to prevent running if data not formatted correctly
  if (NrDims>0 && nrow(RawData)==length(Outcomes)) {

    #construct input vector for search algorithm
    ewInputData = c(NrDims, Iterations, Quantiles)
    for (i in 1:NrDims) {
      ewInputData = c(ewInputData, RawData[i,1:3])
    }
    ewInputData <- as.numeric(ewInputData)

    # set up starting points, upper and lower limits for Z parameters for search algorithm (four numbers are a,b,A,B of a scaled beta where a, b = alpha & beta, A, B = lower and upper limits)
    # hardcoded ranges are (0 - 10, 0 - 10, 0 - 0.99999, 0.00001 - 1) if modCM
    # and (0 - 10, 0 - 10, 0 - 0, 1 - 1) if CM (ie beta dbn, not scaled beta)
    ewStart <- NULL
    ewLower <- NULL
    ewUpper <- NULL
    for (i in 1:NrDims) {
      ewStart <- c(ewStart,1,1,0,1)
      ifelse (ModCMorCM == 1, ewLower <- c(ewLower,0,0,0,0.00001), ewLower <- c(ewLower,0,0,0,1))
      ifelse (ModCMorCM == 1, ewUpper <- c(ewUpper,10, 10,0.99999,1), ewUpper <- c(ewUpper,10,10,0,1))
    }

    #format for crs2lm function:
    # (start point, fn to be minimised, lower bounds, upper bounds, max iterations)
    X <- nloptr::crs2lm(x0 = c(ewInputData,ewStart), TestFitModCM, lower = c(ewInputData,ewLower), upper = c(ewInputData,ewUpper), nl.info=TRUE, xtol_rel=1e-4, maxeval=MaxSearches)
    x <- X$par

    #Pull out a matrix of the Zed parameters from the results
    ZedParams <- matrix(x[(1+length(x)-4*NrDims):(length(x))], nrow = NrDims, ncol = 4, byrow=TRUE, dimnames = list(rep("Z",NrDims),c("a","b","L","U")))

    #make sure scaled betas are properly defined (i.e. L<U)
    #t transposes the output of the apply function (matrix is wrong way around after apply())
    ZedParams <- t(apply(ZedParams, 1, consistentBetas))

    #relabel the columns if they've been swapped using consistentbetas()
    colnames(ZedParams)<- c("a","b","L","U")
    
    #Remove the last row of ZedParams as these are irrelevant (because Z[k]==1)
    ZedParams <- ZedParams[-nrow(ZedParams),]

    #sample from the modCM distribution & find the quantiles
    p <- rModCM(Iterations, ZedParams)
    ewResults <- matrix(quantile(p[,1], Quantiles), nrow = NrDims, ncol = 3, byrow=TRUE)
    for (i in 2:NrDims) {
      ewResults[i,] <- quantile(p[,i], Quantiles)
    }

    #create matrix of results (Tgt=Target from elicitation, Mdl = modelled from ModCM dbn)
    ewSummary <- matrix(rep(0, 6*NrDims), nrow = NrDims, ncol = 6, byrow=TRUE, dimnames = list(rep("P",NrDims),c("Tgt_LL","Tgt_MED","Tgt_UL","Mdl_LL","Mdl_MED","Mdl_UL")))
    for (i in 1:NrDims) {
      ewSummary[i,1:3] <- x[(3+3*i):(5+3*i)]
    }

    ewSummary[1:NrDims,4:6] <- round(ewResults,2)

    #Report coefficients of Zeds, target and actual quantiles
    # Note the coefficients of the last Zed are replaced with Z==1
    ewSummary <- cbind(rbind(ZedParams,c(0,0,0,0)),ewSummary)

    rownames(ewSummary) <- Outcomes

  }
  #print(ewSummary)
  return(ewSummary)
}

#' Fit modified Connor-Mosimann distributions to multiple experts' opinions
#'
#' This function fits a modified CM distribution to an elicited set of quantiles from multiple experts.
#' It saves time by calculating the results for several experts at once, rather than having to call fitModCM() for each expert.
#' It takes a raw datafile with the elicited quantiles from each expert, divides it up expert by expert, and calls
#' fitModCM() once for each expert, returning a data frame with Zed parameters of the modCM distribution for each, along with the target and modelled quantiles.
#'
#' @param Quantiles Sets the quantiles to be fit.  If median and 95\% Credibility Intervals, then set to c(0.025,0.5,0.975).  If median and tertiles then c(0.33,0.5,0.66).  If median and quartiles then c(0.25,0.5,0.75) and so on.
#' @param SearchParams A vector of number of iterations and max number of searches.
#' Number of iterations is the number of draws from the mCM distribution used to estimate the quantiles.  10,000 is a good number to try first.
#' Max number of searches is the maximum number of searches the search algorithm conducts.  The higher the better the fit, but also the longer it takes to find it.  
#' @param RawData Data frame of data elicited from experts.  Must have five columns: expert, outcome, lower quantile, median, upper quantile.
#' @keywords modified Connor-Mosimann
#' @export
#' @examples
#' Quantiles <- c(0.025,0.5,0.975) # to fit median and 95% Credibility Intervals
#' SearchParams <- c(10000,100) # number of iterations, max number of searches
#' RawData <- data.frame(expert = as.character(c(1,1,1,2,2,2)),
#'                       Outcome = as.factor(c("Remission","Progression","Dead",
#'                                            "Remission","Progression","Dead")),
#'                       LL = as.numeric(c(0.43, 0.16, 0.03, 0.35, 0.15, 0.00)),
#'                      MED = as.numeric(c(0.55, 0.27, 0.18, 0.60, 0.30, 0.10)),
#'                       UL = as.numeric(c(0.65, 0.46, 0.23, 0.70, 0.45, 0.20)))
#'
#' fitMultiplemCM(Quantiles, SearchParams, RawData)
#' @return Returns a matrix, each row representing one of the Outcomes for each expert.  Columns a, b, LL and UL are the parameters of the mCM distribution.  
#' Note that the final row will always be zeroes as the number of rows of the mCM distribution is k-1 where k is the number of outcomes.
#' The next three columns (Tgt_LL, Tgt_MED, Tgt_UL) are the target quantiles input to the function as 'RawData' in the example above.
#' The final three (Mdl_LL, Mdl_MED, Mdl_UL) are the quantiles resulting from the model fitting.  If the mCM is a good fit, the Mdl columns will be 
#' identical to the Tgt columns. 
fitMultiplemCM <- function(Quantiles, SearchParams, RawData) {
  Experts <- unique(RawData[,1])

  y <- NULL
  i<- 0
  for (Expert in Experts) {
    i=i+1
    cat("Loop number", i)
    cat(Expert)
    ewData <- subset(RawData,RawData[,1]==Expert)
    x <- fitModCM(ewData[,2], ewData[,3:5], SearchParams, 1, Quantiles)
    x <- cbind(expert=Expert,x)
    y <-rbind(y,x)
  }

  #tidy up y and return as a data frame
  y <- cbind(y[,1],rownames(y),y[,2:ncol(y)])
  rownames(y) <- NULL
  colnames(y)[1:6] <- c("Expert","Outcome","a","b","L","U")
  y <- as.data.frame(y)
  y[,3:ncol(y)] <- apply(y[,3:ncol(y)],2,as.numeric)
  return(y)
}

#' Merge multiple modified Connor-Mosimann distributions together
#'
#' This function merges the results of multiple experts' distributions using a numeric linear pool approach.
#' It samples from the distributions of all experts individually many times (e.g. 100,000), then calculates the overall quantiles and medians from the samples.
#' The function returns a matrix representing the lower, median and upper limits of the pooled distribution.
#' This can then be fed into fitModCM() to generate a modified Connor-Mosimann distribution representing the overall spread of the experts' beliefs.
#'
#' @param NrSamples Vector of length 1.  Sets the number of samples to draw from each expert's mCM distribution
#' @param RawData Data frame of mCM parameters.  Must have six columns: expert, outcome, a, b, L, U of modified Connor-Mosimann distribution.  Note last row of parameters will always be zero.
#' Columns 1:6 of the output from function fitMultipleCM() are in the correct format for this.
#' @keywords modified Connor-Mosimann
#' @export
#' @examples
#' NrSamples <- 100000
#' RawData <- data.frame(expert = as.character(c(1,1,1,2,2,2)),
#'                       Outcome = as.factor(c("Remission","Progression","Dead",
#'                                            "Remission","Progression","Dead")),
#'                       a = as.numeric(c(6.0786, 0.2245, 0, 6.9214, 4.5259, 0)),
#'                       b = as.numeric(c(7.5900, 0.5866, 0, 1.7187, 3.1892, 0)),
#'                       L = as.numeric(c(0.3400, 0.4839, 0, 0.0152, 0.2390, 0)),
#'                       U = as.numeric(c(0.7917, 0.9213, 0, 0.7106, 0.9970, 0)))
#' mergeMultiplemCM(NrSamples,RawData)
mergeMultiplemCM <- function(NrSamples,RawData) {
  #first create comprehensive list of all experts and outcomes
  Experts <- unique(RawData[,1])
  outcomes <- unique(RawData[,2])

  #create matrix to hold overall results
  pAll <- matrix(data=0,nrow=0,ncol=length(outcomes))
  colnames(pAll) <- outcomes

  # Now sample from each distribution
  for (Expert in Experts) {
    cat("\n Expert",Expert)
    #read in the dataExpert
    x <- RawData[RawData[,1]==Expert,]

    #sample from it
    p <- rModCM(NrSamples,x[-nrow(x),3:6])
    colnames(p)<- x[,2]
    #add results to combined results
    #first lengthen pAll by NrSamples
    pAll <- rbind(pAll,matrix(data=0,nrow=NrSamples,ncol=ncol(pAll)))

    #now work through each column in p and match it to equivalent in pAll
    for (j in 1:ncol(p)) {
      pAll[(nrow(pAll)-NrSamples+1):nrow(pAll),which(colnames(pAll)==colnames(p)[j])] = p[,j]
    }
  }

#calculate the median and upper and lower CrIs
  ewQuantiles <- matrix(quantile(pAll[,1], c(0.025,0.5,0.975)), nrow = ncol(pAll), ncol = 3, byrow=TRUE)
  for (i in 2:ncol(pAll)) {
    ewQuantiles[i,] <- quantile(pAll[,i], c(0.025,0.5,0.975))
  }
  rownames(ewQuantiles) <- outcomes
  colnames(ewQuantiles) <- c("LL","MED","UL")

#sort in descending order of medians
  ewQuantiles <- ewQuantiles[order(ewQuantiles[,2],decreasing=TRUE),]
  ewQuantiles <- as.matrix(ewQuantiles)
  return(ewQuantiles)
}



#' Fit Dirichlet distribution
#'
#' This function fits a Dirichlet distribution to an elicited set of quantiles from, e.g. an expert elicitation workshop.
#' The function uses the crs2lm() function from the nloptr package to search for the
#' set of hyperparameters that that generates quantiles that match the elicited data as closely as possible.
#' crs2lm() repeatedly calls TestFitDirichlet() until it finds the best fitting set of inputs.
#' The function returns a data frame containing the best fit set of hyperparameters of the Dirichlet, the Target quantiles and the modelled quantiles.
#'
#' @param Outcomes A vector of names of outcomes
#' @param RawData matrix of lower, middle and upper quantiles for each dimension elicited from experts
#' @param SearchParams A vector of number of iterations and max number of searches.
#' Number of iterations is the number of draws from the mCM distribution used to estimate the quantiles.  Try 10,000 first.
#' Max number of searches is the maximum number of searches the search algorithm conducts.  The higher the better the solution, but also the longer it takes.  Try 1000 first.
#' @param Quantiles Sets the quantiles to be fit.  If median and 95\% Credibility Intervals, then set to c(0.025,0.5,0.975).  If median and tertiles then c(0.33,0.5,0.66).  If median and quartiles then c(0.25,0.5,0.75) and so on.
#' @keywords Dirichlet
#' @export
#' @examples
#' Outcomes <- c("Remission","Progression","Dead")
#' RawData <- matrix(data = c(0.43, 0.55, 0.65,
#'                           0.16, 0.27, 0.46,
#'                           0.03, 0.18, 0.23
#'            ),ncol=3,byrow=TRUE)
#'
#' SearchParams <- c(10000,100) #number of iterations, max number of searches
#'
#' Quantiles <- c(0.025,0.5,0.975) # example here is 95% credibility limits and median.
#'
#' fitDirichlet(Outcomes, RawData, SearchParams, Quantiles)
#' @return Returns a matrix, each row representing one of the Outcomes.  Column 'Dirichlet' is the parameters of the Dirichlet distribution.  
#' The next three columns (Tgt_LL, Tgt_MED, Tgt_UL) are the target quantiles input to the function as 'RawData' in the example above.
#' The final three (Mdl_LL, Mdl_MED, Mdl_UL) are the quantiles resulting from the model fitting.  If the mCM is a good fit, the Mdl columns will be 
#' identical to the Tgt columns. 
fitDirichlet <- function(Outcomes, RawData, SearchParams, Quantiles) {
  
  NrDims <- length(Outcomes)
  Iterations <- SearchParams[1]
  MaxSearches <- SearchParams[2]
  
  #If statement to prevent running if data not formatted correctly
  if (NrDims>0 && nrow(RawData)==length(Outcomes)) {
    
    #Format for input data to search algorithm:
    # number of dimensions to distribution; FromState; ToState1; LL1; MED1; UL1; ToState2; LL2; MED2; UL2... etc
    ewInputData = c(NrDims, Iterations, Quantiles)
    for (i in 1:NrDims) {
      ewInputData = c(ewInputData, RawData[i,1:3])
    }
    ewInputData <- as.numeric(ewInputData)
    
    # set up starting points, upper and lower limits for Dirichlet parameters for search algorithm
    ewStart <- rep(1,NrDims)
    ewLower <- rep(0,NrDims)
    ewUpper <- rep(1000000,NrDims)
    
    #testlines
    #print("Inputs:")
    #print(ewStart)
    #print(ewLower)
    #print(ewUpper)
    #end testlines
    
    #format for crs2lm function:
    # (start point, fn to be minimised, lower bounds, upper bounds, max iterations)
    X <- nloptr::crs2lm(x0 = c(ewInputData,ewStart), TestFitDirichlet, lower = c(ewInputData,ewLower), upper = c(ewInputData,ewUpper), nl.info=TRUE, xtol_rel=1e-4, maxeval=MaxSearches)
    x <- X$par
    
    #Pull out a matrix of the dirichlet parameters from the results
    DirichletParams <- matrix(x[(length(x)-NrDims+1):(length(x))], nrow = NrDims, byrow=TRUE)
    colnames(DirichletParams)<-"Dirichlet"
    
    p <- gtools::rdirichlet(Iterations,DirichletParams)
    
    #Find the quantiles from the sampled values
    ewResults <- matrix(quantile(p[,1], Quantiles), nrow = NrDims, ncol = 3, byrow=TRUE)
    for (i in 2:NrDims) {
      ewResults[i,] <- quantile(p[,i], Quantiles)
    }
    
    #create matrix of results (Tgt=Target from elicitation, Mdl = modelled from ModCM dbn)
    
    ewSummary <- matrix(rep(0, 6*NrDims), nrow = NrDims, ncol = 6, byrow=TRUE, dimnames = list(rep("P",NrDims),c("Tgt_LL","Tgt_MED","Tgt_UL","Mdl_LL","Mdl_MED","Mdl_UL")))
    for (i in 1:NrDims) {
      ewSummary[i,1:3] <- x[(3+3*i):(5+3*i)]
    }
    
    ewSummary[1:NrDims,4:6] <- round(ewResults,4)
    
    #Report coefficients of Zeds, target and actual quantiles
    # Note the coefficients of the last Zed are replaced with Z==1
    ewSummary <- cbind(DirichletParams,ewSummary)
    
    rownames(ewSummary) <- Outcomes
    
  }
  print(ewSummary)
  return(ewSummary)
}


# Calculate sum of squared deviation between modelled and elicited quantiles for Dirichlet distribution
#
# This function uses a set of hyperparameters of a Dirichlet distribution as an input and calculates the
# associated quantiles of the marginal betas for that distribution, then compares this
# with the target quantiles by calculating the sum of squared deviation.
# The function returns the SSD.
# @param x Vector of input parameters: NrDims, Iterations, Quantiles being compared,
# Target for quantiles and hyperparameters of Dirichlet distribution.
# Note input must be a single atomic vector not dataframes or matrices.  This is because the
# search algorithm crs2lm(), called in fitDirichlet() can only send a single vector when it calls this function.
# @param NrDims Number of dimensions to mCM distribution (e.g. number of outcomes such as remission, progression, dead)
# @param Iterations Number of Samples from the mCM distribution to pull to calculate the quantiles.  The higher the more precise the reported quantiles.
# @param Quantiles Vector of quantiles to be matching.
# @param Targets Target values for lower, middle and upper quantiles for each dimension (i.e. quantities elicited from experts).  Order of numbers is LL1, Median1, UL1, LL2, Median2, UL2 etc
# @param Dirich Hyperparamters of Dirichlet distribution from which to calculate the quantiles.
# @keywords Dirichlet
# @examples
# NrDims <- 3
# Iterations <-10000
# Quantiles <- c(0.025,0.5,0.975) #(median and upper and lower 95% credibility interval)
# Targets <- c(0.43, 0.55, 0.65,
#              0.16, 0.27, 0.46,
#              0.03, 0.18, 0.23)
# Dirich <- c(10,10,10)
# TestFitDirichlet(c(NrDims,Iterations,Quantiles,Targets,Dirich))
#
# 
TestFitDirichlet <- function(x) {
  
  x <- as.numeric(x)
  
  #testline
  #print(x)
  #end testline
  
    #format of x:
  #NrDims, Iterations, Quantiles
  #L1,M1,U1,L2,M2,U2,L3,M3,U3,...
  #n1, n2, n3....
  
  #assemble input vector into something easier to understand
  NrDims <- x[1]
  Iterations <- x[2]
  Quantiles <- x[3:5]
  
  #target LL, MED, ULs for each dimension
  ewTargets <- matrix(data=x[6:(5+3*NrDims)], nrow = NrDims, ncol = 3, byrow = TRUE)
  
  #Dirichlet parameter values
  Dirich <- x[(6+3*NrDims):length(x)]
  
  #testlines
  #print(paste("NrDims",NrDims))
  #print(paste("Iterations", Iterations))
  #print("Quantiles")
  #print(Quantiles)
  #print("Targets")
  #print(ewTargets)
  #print("Dirich")
  #print(Dirich)
  #print("sum(Dirich)")
  #print(sum(Dirich))
  #end testlines
  
  #error trap in case any Dirich parameters specified as <0.1
  if (sum(Dirich<0.1)>0) Dirich <- rep(1,length(Dirich))

  #testline
  #print(Dirich)
  #end testlines
    
  #calculate the p's
  p <- gtools::rdirichlet(Iterations,Dirich)
  
  #Find the quantiles from the sampled values
  ewResults <- matrix(quantile(p[,1], Quantiles), nrow = NrDims, ncol = 3, byrow=TRUE)
  for (i in 2:NrDims) {
    ewResults[i,] <- quantile(p[,i], Quantiles)
  }
  
  #Code to place greater weight on getting the medians right over the LL & ULs
  #ewWeights <- matrix(rep(c(1,NrDims,1),NrDims), nrow = NrDims, ncol = 3, byrow=TRUE)
  #unweighted
  ewWeights <- matrix(rep(c(1,1,1),NrDims), nrow = NrDims, ncol = 3, byrow=TRUE)
  
  
  #calculate (weighted) sum of the squared deviation from target
  SSD = sum((ewWeights*(ewTargets-ewResults))^2)
  cat("\n",SSD)
  return(SSD)
}
