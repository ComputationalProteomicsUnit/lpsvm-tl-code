#-----------------------------------------------------------------------------
# Linear programming SVM as described in Wu and Dietterich, "Improving SVM 
# accuracy by Training on Auxilliary Data Sources".
#
# The only difference is that this enforces xi as the unconstrained 
# version does not appear to work.
#
# This version is modified to work with two kernels, generally one for
# each of a pair of data sources with different features.
#
# Copyright (C) Sean Holden 2014-15.
#-----------------------------------------------------------------------------

source("kernels-c.R")
require("lpSolve")

#-----------------------------------------------------------------------------
# Learn from the data by setting up and solving the linear program. This 
# version implements the simplest version of adding a second kernel 
# for auxilliary data, but should be useable in the normal case exactly as 
# with the earlier version.
#
# data contains X and y as in the standard version. C, gamma and balanced are
# also the same as for the standard version. To use two kernels set twoK
# to TRUE. Then X2 contains the extra vectors for the second kernel in row
# to row correspondence with data$X, and gamma2 is the parameter for the \
# second RBF kernel.
#-----------------------------------------------------------------------------
lpSVMLearn <- function(data, C, gamma, balanced = FALSE, twoK = FALSE, 
              X2 = matrix(), gamma2 = 1) {
 
 	X <- data$X
 	y <- data$y
 	stopifnot((1 %in% y) && (-1 %in% y))
    m = nrow(X)
    d = ncol(X)
   
    cVector = rep(C, times = m)
    # Deal with an unbalanced data set.
    if (!balanced) {
        positives <- (y == 1)
        negatives <- (y == -1)
        nplus <- sum(positives)
        nminus <- m - nplus
        cVector[positives] <- cVector[positives] * sqrt(nminus/nplus)
        cVector[negatives] <- cVector[negatives] * sqrt(nplus/nminus)
    }
    
    # Need to construct the linear program.
	# First the constant vector defining what's to be minimized.
    # The variable vector has alpha, b, xi (assuming 1 kernel). 
    c <- c(rep(1, times = m), rep(0, times = 2), cVector)

   	# Now the tricky bit. We need to make the constraint matrix and vector.
   	# The vector is straightforward.
   	b <- rep(-1, times = m)

   	# The matrix needs several stages.
   	outer <- (y %o% y)
   	B <- gram(X, gamma) * outer

    # And we now check to see if there's a second kernel.
    if (!twoK) 
		A <- cbind(-B, -y, y, -diag(m))
	else {
		c <- c(rep(1, times = m), c)
		B2 <- gram(X2, gamma2) * outer 
		A <- cbind(-B, -B2, -y, y, -diag(m))
	}
	
	# Solve the underlying optimization problem.
    result <- lp("min", c, A, rep("<=", nrow(A)), b)

    # Reconstruct the correct alpha values. 
    r <- result$solution
	# b needs to be negated as Mangasarian uses its negative.
	if (!twoK) {
		return(list(alpha = r[1:m], b = -(r[m+1] - r[m+2]), 
				    status = result$status, value = result$objval))
    }
    else {
    	twom <- 2 * m
    	return(list(alpha = r[1:m], alpha2 = r[(m+1):twom], 
    				b = -(r[twom+1] - r[twom+2]), 
    		        status = result$status, value = result$objval))
    }
}

#-------------------------------------------------------------------------------
# Given the result of lpSVMLearn and the kernel parameter, compute
# the output of the SVM as a real for a new input x.
#
# Now if twoK is TRUE, X2 and gamma2 are as for lpSVMLearn, and x2 is the
# part of the new x to be classified corresponding to X2. (The actual vector
# you are classifying is the concatenation of x and x2.)
#-------------------------------------------------------------------------------
lpSVM <- function(parameters, data, gamma, x, twoK = FALSE, 
              X2 = matrix(), gamma2 = 1, x2 = c()) {
  	Xt <- t(data$X)
  	m <- ncol(Xt)
	y <- data$y  
    k <- .C("makeVectorRBF", as.double(Xt), as.double(x), 
			as.integer(m), as.integer(nrow(Xt)), as.double(gamma), 
			result = rep(0, m))
    single <- (t(k$result) %*% (y * parameters$alpha)) - parameters$b
    if (!twoK) return(single)
    else {
    	X2t <- t(X2)
    	m2 <- ncol(X2t)
    	k2 <- .C("makeVectorRBF", as.double(X2t), as.double(x2), 
			as.integer(m2), as.integer(nrow(X2t)), as.double(gamma2), 
			result = rep(0, m2))
        return((t(k2$result) %*% (y * parameters$alpha2)) + single)
    }
        
}

#-------------------------------------------------------------------------------
# Same as lpSVM but returns a function computing the output for a given x (and
# if necessary x2).
#-------------------------------------------------------------------------------
lpSVMf <- function(parameters, data, gamma, twoK = FALSE, 
              X2 = matrix(), gamma2 = 1) {
              	
	if (!twoK) f <- function(x) { lpSVM(parameters, data, gamma, x) }
	else {
		f <- function(x, x2) {
			        lpSVM(parameters, data, gamma, x, twoK, X2, gamma2, x2)
			        }
		}
	return(f)
}

