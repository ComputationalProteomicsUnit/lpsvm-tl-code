#------------------------------------------------------------------------
# Kernels for SVMs, and a function to compute the Gram matrix.
#
# Copyright (C) Sean Holden 2014-15.
#------------------------------------------------------------------------

dyn.load("c/outer.so")

#------------------------------------------------------------------------
# Make a Gram matrix from data and kernel.
# The kernel computed is the standard RBF with the supplied parameter
# gamma. X should have one example per row.
#------------------------------------------------------------------------
gram <- function(X, gamma) {
	Xt <- t(X)
	m <- nrow(X)
	g <- .C("makeGramRBF", as.double(Xt), as.integer(m), 
		 as.integer(ncol(X)), as.double(gamma),
		 result = matrix(0,m,m))
	g$result
}


