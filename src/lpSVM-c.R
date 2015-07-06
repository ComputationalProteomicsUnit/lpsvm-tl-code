#-----------------------------------------------------------------------------
# Linear programming SVM as described in Wu and Dietterich, "Improving SVM 
# accuracy by Training on Auxilliary Data Sources".
#
# The only difference is that this enforces xi as the unconstrained 
# version does not appear to work.
#
# Copyright (C) Sean Holden 2014-15.
#-----------------------------------------------------------------------------

source("kernels-c.R")
require("lpSolve")

#-----------------------------------------------------------------------------
# Learn from the data by setting up and solving the linear program.
#
# The kernel K is always RBF at present, and you need to supply gamma.
# Data is assumed to be unbalanced by default.
# C is the usual parameter.
#-----------------------------------------------------------------------------
lpSVMLearn <- function(data, C, gamma, balanced = FALSE) {
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
    # The variable vector has alpha, b, xi. 
    c <- c(rep(1, times = m), rep(0, times = 2), cVector)

    # Now the tricky bit. We need to make the constraint matrix and vector.
    # The vector is straightforward.
    b <- rep(-1, times = m)

    # The matrix needs several stages.
    B <- gram(X, gamma) * (y %o% y)
    A <- cbind(-B, -y, y, -diag(m))

    # Solve the underlying optimization problem.
    result <- lp("min", c, A, rep("<=", nrow(A)), b)

    # Reconstruct the correct alpha values. 
    r <- result$solution
    
    # b needs to be negated as Mangasarian uses its negative.
    list(alpha = r[1:m], b = -(r[m+1] - r[m+2]), status = result$status,
         value = result$objval)
}

#-------------------------------------------------------------------------------
# Given the result of lpSVMLearn and the gamma parameter, compute
# the output of the SVM as a real for a new input x.
#-------------------------------------------------------------------------------
lpSVM <- function(parameters, data, gamma, x) {
	Xt <- t(data$X)
	m <- ncol(Xt)
	k <- .C("makeVectorRBF", as.double(Xt), as.double(x), 
			as.integer(m), as.integer(nrow(Xt)), as.double(gamma), 
			result = rep(0, m))
    (t(k$result) %*% (data$y * parameters$alpha)) - parameters$b    
}

#-------------------------------------------------------------------------------
# Same as lpSVM but returns a function computing the output for a given x.
#-------------------------------------------------------------------------------
lpSVMf <- function(parameters, data, gamma) {
	function(x) { lpSVM(parameters, data, gamma, x) }
}

