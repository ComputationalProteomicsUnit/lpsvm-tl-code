#-------------------------------------------------------------------------------
# Linear programming SVM as described in Mangasarian, Generalized Support
# Vector Machines, section 4 equation (11) and with probabilistic output
# as described by Chang and Lin, LIBSVM: A Library for Support Vector 
# Machines.
#
# This version is modified to work with two kernels, generally one for
# each of a pair of data sources with different features.
#
# Copyright (C) Sean Holden 2014-15.
#-------------------------------------------------------------------------------

source("lpSVM2K-c.R")
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

#-------------------------------------------------------------------------------
# As lpSVMLearn but using probabilistic outputs.
#
# fs is the number of folds used for probability estimation.
#
# For a description of the other parameters see lpSVMlearn in the file
# lpSVM2K-c.R
#-------------------------------------------------------------------------------
lpSVMLearnP <- function(data, C, gamma, balanced = FALSE,  
              X2, gamma2, fs = 5) {
	
	X <- data$X
	y <- data$y
	m <- nrow(X)
	d <- ncol(X)
	d2 <- ncol(X2)
	
	stopifnot((1 %in% y) && (-1 %in% y))
	
	# To get the probability coefficients we need to get some 
	# cross-validated outputs. This means joining the data together.
	# I'm now going to drop the option of keeping this compatible with 
	# the single kernel case as it's getting horribly tricky.
	joinedData <- cbind(X, X2)
	folds <- makeStratifiedCVFolds(list(X=joinedData,y=y), fs)
	cvf <- NULL
	cvy <- NULL
	for (i in 1:fs) {
		currentSets <- cvMakeSets(folds, i)
		thisX <- currentSets$X
		thisy <- currentSets$y
		thisTestX <- currentSets$testX
		parameters <- lpSVMLearn(list(X=thisX[,1:d], y=thisy), C, gamma, 
								 balanced, TRUE, thisX[,(d+1):(d+d2)], gamma2)
		app <- c()
		for (ex in 1:nrow(thisTestX)) 
			app[ex] <- lpSVM(parameters, list(X = thisX[,1:d], y = thisy), 
						gamma, thisTestX[ex, 1:d], TRUE, thisX[,(d+1):(d+d2)], 
				 		gamma2, thisTestX[ex, (d+1):(d+d2)])
		cvf <- c(cvf, app)
		cvy <- c(cvy, currentSets$testy)
	}
	coeffs <- estimateProbabilityCoeffs(cvf, cvy)
	
	# Finally, do the learning and use the cross-validated values to find the 
	# parameters to obtain probabilities.
	parameters <- lpSVMLearn(data, C, gamma, balanced, TRUE, X2, gamma2)
	parameters$pcoeffs <- coeffs
	return(parameters)
}

#-------------------------------------------------------------------------------
# As lpSVM, but with probabilistic output.
#-------------------------------------------------------------------------------
lpSVMP <- function(parameters, data, gamma, x, X2, gamma2, x2) {
	1/(1 + exp(((parameters$pcoeffs$A) * 
	lpSVM(parameters, data, gamma, x, TRUE, X2, gamma2, x2)) + 
	parameters$pcoeffs$B))
}

#-------------------------------------------------------------------------------
# As lpSVMf, but with probabilistic output.
#-------------------------------------------------------------------------------
lpSVMPf <- function(parameters, data, gamma, X2, gamma2) {
	function(x, x2) { lpSVMP(parameters, data, gamma, x, X2, gamma2, x2) }
}

