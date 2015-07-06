#-----------------------------------------------------------------------------
# Linear programming SVM as described in Wu and Dietterich, "Improving SVM 
# accuracy by Training on Auxilliary Data Sources".
#
# The only difference is that (a) this enforces xi as the unconstrained 
# version does not appear to work and (b) this implements multiclass 
# classification in the same way as libSVM.
#
# This version is modified to work with two kernels, generally one for
# each of a pair of data sources with different features.
#
# Copyright (C) Sean Holden 2014-15.
#-----------------------------------------------------------------------------

source("lpSVM2K-c.R")
source("multiShared.R")

#-------------------------------------------------------------------------------
# Linear programming SVM learning for multiple-class 
# problems. Uses the binary version lpSVMLearn.
#
# y should specify classes as positive integers.
#
# For a description of the parameters see lpSVMlearn in the file
# lpSVM2K-c.R
#-------------------------------------------------------------------------------
lpSVMLearnMulti <- function(data, X2, C, gamma, gamma2, balanced = FALSE) {
	y <- data$y
	stopifnot(length(unique(y)) > 1)
	
	# Need to join the X and X2 data together to split into binary problems.
	X <- data$X
	c1 <- ncol(X)
	c2 <- ncol(X2)
	binary <- multiToBinary(list(X = cbind(X, X2), y = y))
	numBinary <- length(binary)
	
	svmLearn <- function(d, x2) { lpSVMLearn(d, C, gamma, balanced, 
									twoK = TRUE, x2, gamma2) }
	
	result <- list()
	for (i in 1:numBinary) {
		# Need to split the X data again.
		binAllX <- binary[[i]]$X
		binX2 <- binAllX[,(c1+1):(c1+c2)]
		binary[[i]]$X <- binAllX[,1:c1]  
		binary[[i]]$X2 <- binX2
		currentParameters <- svmLearn(binary[[i]], binX2)
		# The result needs to maintain the training set used as this is
		# needed for classification.
		result[[i]] <- list(parameters = currentParameters, 
							data=binary[[i]])
    	}
	return(result)
}

#-------------------------------------------------------------------------------
# Given the parameters from training a multiclass lpSVM with lpSVMLearnMulti
# classify a new example x, x2.
#-------------------------------------------------------------------------------
lpSVMMulti <- function(pList, y, gamma, gamma2, x, x2) {
	numBinary <- length(pList)
	classes <- unique(y)
	maxClass <- max(classes)
	
	# First make a vector of individual classifications.
	# The parametersList contains for each binary problem a pair denoting 
	# which two classes it addresses.
	output <- rep(0, numBinary)	
	for (i in 1:numBinary) {
		p <- pList[[i]]

		# If the linear program was solvable for the current 
		# problem then find the predicted class.
		if (p$parameters$status == 0) {
			out <- lpSVM(p$parameters, p$data, gamma, x, twoK = TRUE, p$data$X2, 
					 gamma2, x2) 
			if (out >= 0) output[i] <- p$data$classes[[1]] 
			         else output[i] <- p$data$classes[[2]]	
		}
		# If the linear program was *not* solvable then predict the 
		# most abundant class from the training data.
		else {
			classOne <- sum(p$data$y == 1)
			classTwo <- sum(p$data$y == -1)
			if (classOne >= classTwo) 
				output[i] <- p$data$classes[[1]]
			else
				output[i] <- p$data$classes[[2]]
		}	
	}
	
    # Now see which class has won.
   	counts <- rep(0, maxClass)
	for (i in classes) counts[i] <- sum(output == i)
	winner <- max(counts)
	fullOutput<-(1:maxClass)[counts == winner]
	return(fullOutput[1])          # Which is dumb but what libSVM does.
}

#-------------------------------------------------------------------------------
# Return a function for classifying new examples.
#
# d1 is the number of columns covering the first chunk of X data.
#-------------------------------------------------------------------------------
lpSVMMultif <- function(plist, d1, y, gamma, gamma2) {
	function(input) { 
		x <- input[1:d1]
		x2 <- input[(d1+1):length(input)]
		lpSVMMulti(plist, y, gamma, gamma2, x, x2) 
	}
} 
