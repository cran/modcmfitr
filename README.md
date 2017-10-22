modcmfitr
=========

The purpose of this package is to fit a modified Connor-Mosimann distribution to quantiles for multinomial problems elicited from experts. It also includes functions to fit Connor-Mosimann and Dirichlet distributions, and to sample from modified Connor-Mosimann and Connor-Mosimann distributions for use in decision analytic models (rtools already provides a function to sample from a Dirichlet distribution).

For details on how to use please see the accompanying vignette.

Example
-------

``` r

# Fit a distribution to elicited quantiles
Outcomes <- c("Remission","Progression","Dead")
RawData <- matrix(data = c(0.43, 0.55, 0.65,
                         0.16, 0.27, 0.46,
                         0.03, 0.18, 0.23
          ),ncol=3,byrow=TRUE)
SearchParams <- c(10000,100) #number of iterations, max number of searches
ModCMorCM <- 1 # if 1 will fit mCM, if 0 will fit CM
Quantiles <- c(0.025,0.5,0.975) # example here is 95% credibility limits and median.
mCM <- fitModCM(Outcomes, RawData, SearchParams, ModCMorCM, Quantiles)

# Sample from the fitted distribution
n <- 100
Z <- mCM[1:(nrow(mCM)-1),1:4]
mCMSamples <- rModCM(n,Z)
colnames(mCMSamples) <- Outcomes
```

Installation instructions
-------------------------

``` r
  install.packages("modcmfitr")
```
