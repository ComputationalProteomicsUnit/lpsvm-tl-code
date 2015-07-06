## Linear Programming SVM Transfer Learning
## https://github.com/ComputationalProteomicsUnit/lpsvm-tl-code
## Copyright (C) 2014-2015  Sean Holden
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License along
## with this program; if not, write to the Free Software Foundation, Inc.,
## 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

#-------------------------------------------------------------------------------
# Linear programming SVM as described in Mangasarian, Generalized Support
# Vector Machines, section 4 equation (11) and with probabilistic output
# as described by Chang and Lin, LIBSVM: A Library for Support Vector 
# Machines.
#
#-------------------------------------------------------------------------------

source("lpSVM-c.R")
source("cv.R")

#-----------------------------------------------------------------------------
# Estimate the coefficients for interpreting the output of an SVM
# probabilistically using the algorithm from Lin et al, A Note on
# Platt's Probabilistic Outputs for Support Vector Machines. Directly
# follows the pseudocode given in Appendix 3.
#
# f is the SVM decision values
# y is the training labels
#-----------------------------------------------------------------------------
estimateProbabilityCoeffs <- function(f, y, maxIter = 100, verbose = FALSE,
                                minStep = 1e-10, sigma = 1e-12) {
    failed <- FALSE 
    maxIterReached = FALSE  	
       								  	
    labels <- y
    labels[y == -1] <- 0   	
    prior1 <- sum(y == 1)
    prior0 <- sum(y == -1)	
    
    # Initial values
    hiTarget <- (prior1 + 1)/(prior1 + 2)
    loTarget <- 1/(prior0 + 2)
    l <- prior1 + prior0
    t <- rep(loTarget, l)
    t[labels > 0] <- hiTarget
	
    a <- 0
    b <- log((prior0 + 1)/(prior1 + 1))
    fval <- 0	
    for (i in 1:l) {
      fApB <- ((f[i] * a) + b)
      if (fApB >= 0)
        fval <- fval + (t[i] * fApB) + log(1 + exp(-fApB))
      else
        fval <- fval + ((t[i] - 1) * fApB) + log(1 + exp(fApB))
    }		
    
    for (it in 1:maxIter) {
      if (verbose) cat("Iteration:",it,"\n")
      h11 <- sigma
      h22 <- sigma
      h21 <- 0
      g1 <- 0
      g2 <- 0
      for (i in 1:l) {
        fApB <- f[i]*a + b
        if (fApB >= 0) {
          q <- 1/(1 + exp(-fApB))
    	  p <- exp(-fApB) * q 
    	}
    	else {
    	  p <- 1/(1+exp(fApB))
    	  q <- p * exp(fApB)
    	}
    	d2 <- p * q
    	h11 <- h11 + f[i] * f[i] * d2
    	h22 <- h22 + d2
    	h21 <- h21 + f[i] * d2
    	d1 <- t[i] - p
    	g1 <- g1 + f[i] * d1
    	g2 <- g2 + d1
      }
      if ((abs(g1) < 1e-5) && (abs(g2) < 1e-5)) break
      det <- h11 * h22 - h21 * h21
      dA <- -(h22 * g1 - h21 * g2)/det
      dB <- -(-h21 * g1 + h11 * g2)/det
      gd <- g1 * dA + g2 * dB
      stepSize <- 1
      while (stepSize >= minStep) {
        if (verbose) cat("stepSize =",stepSize,"\n")
    	newA <- a + stepSize * dA
    	newB <- b + stepSize * dB
    	newf <- 0
    	for (i in 1:l) {
    	  fApB <- f[i] * newA + newB
    	  if (fApB >= 0)
    	    newf <- newf + t[i] * fApB + log(1 + exp(-fApB))
    	  else
            newf <- newf + (t[i] - 1) * fApB + log(1 + exp(fApB))
    	}
    	if (newf < fval + 0.0001 * stepSize * gd) {
    	  a <- newA
    	  b <- newB
    	  if (verbose) cat("A =",a,"B =",b,"\n")
    	  fval <- newf
    	  break	
    	}
    	else {
    	  stepSize <- stepSize / 2
    	}
    	if (stepSize < minStep) {
    	  failed <- TRUE
          break
    	}
      }
      if (it >= maxIter) maxIterReached <- TRUE
    }	
    return(list(A=a, B=b, status=failed, max=maxIterReached))	
}

#-----------------------------------------------------------------------------
# As lpSVMLearn but using probabilistic outputs.
#
# fs is the number of cross-validation folds used for probability estimation.
#-----------------------------------------------------------------------------
lpSVMLearnP <- function(data, C, gamma, balanced = FALSE, fs = 5) {
  y <- data$y
  stopifnot((1 %in% y) && (-1 %in% y))
	
  # To get the probability coefficients we need to get some 
  # cross-validated outputs.
  folds <- makeStratifiedCVFolds(data, fs)
  cvf <- NULL
  cvy <- NULL
  for (i in 1:fs) {
    currentSets <- cvMakeSets(folds, i)
    parameters <- lpSVMLearn(currentSets, C, gamma, balanced)
    cvf <- c(cvf, apply(currentSets$testX, 1, 
			lpSVMf(parameters, currentSets, gamma)))
    cvy <- c(cvy, currentSets$testy)
  }
  coeffs <- estimateProbabilityCoeffs(cvf, cvy)
	
  # Finally, do the learning and use the cross-validated values to find the 
  # parameters to obtain probabilities.
  parameters <- lpSVMLearn(data, C, gamma, balanced)
  parameters$pcoeffs <- coeffs
  return(parameters)
}

#-------------------------------------------------------------------------------
# As lpSVM, but with probabilistic output.
#-------------------------------------------------------------------------------
lpSVMP <- function(parameters, data, gamma, x) {
  1/(1 + exp(((parameters$pcoeffs$A) * lpSVM(parameters, data, gamma, x)) + 
     parameters$pcoeffs$B))
}

#-------------------------------------------------------------------------------
# As lpSVMf, but with probabilistic output.
#-------------------------------------------------------------------------------
lpSVMPf <- function(parameters, data, gamma) {
	function(x) { lpSVMP(parameters, data, gamma, x) }
}

