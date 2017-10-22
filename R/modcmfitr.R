#' modcmfitr: A package to fit and sample from a modified Connor-Mosimann distribution
#'
#' This package fits a modified Connor-Mosimann (mCM) distribution to quantiles of a multinomial distribution
#' elicited from experts using a process such as the Sheffield Elicitation Framework (SHELF, O'Hagan & Oakley). 
#' 
#' More details and worked examples are provided in the help files and accompanying vignette.  Functions are also provided to fit Connor-Mosimann (CM) and Dirichlet (D) distributions, and to sample from mCM and CM
#' distributions.  A function to sample from the Dirichlet distribution is available in the gtools package (Warnes, Bolker & Lumley).
#'
#' Copyright (C) 2017  Ed Wilson, University of Cambridge, <ed.wilson@medschl.cam.ac.uk>
#' 
#' October 2017
#' 
#' @docType package
#' @name modcmfitr
#' @importFrom stats quantile
#' @importFrom stats qbeta
#' @importFrom stats rbeta
#' 
NULL