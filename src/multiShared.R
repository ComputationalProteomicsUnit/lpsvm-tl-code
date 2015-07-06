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

#-----------------------------------------------------------------------------
# Utilities employed by the implementation of multiple-class linear
# programming SVMs.
#
#-----------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Given pairs corresponding to attribute vectors Xi and their corresponding 
# label yi, construct a set of binary problems, one for each class pair. 
#
# The output uses the order (1,1), (1,2), ..., (1,n), (2,2), (2,3), where 
# the numbers here represent indexes of the ys parameter. Returned results are 
# labelled with the actual pair of labels they correspond to.
#-------------------------------------------------------------------------------
pairClasses <- function(Xs, ys) {
    n <- length(ys)
    k <- 1
    result <- list()
    for (i in 1:(n-1)) {
    	for (j in (i+1):n) {
    		c1 <- ys[i]
    		c2 <- ys[j]
            result[[k]] <- list(Xs[[c1]], Xs[[c2]], c1, c2)
            k <- k + 1
        }
    }
    return(result)
}

#-------------------------------------------------------------------------------
# From multiclass data, make binary problems with labels converted to
# +1/-1.  All pairs of classes, with corresponding binary training/test sets
# are returned in a list. Labels supplied should be positive integers.
#
# The parameter randomize causes the order to be randomly permuted.
# Alternatively ybin has all +1 followed by all -1.
#-------------------------------------------------------------------------------
multiToBinary <- function(data, randomize=TRUE) {
	X <- data$X
	y <- data$y

    labels <- unique(y)
        
    # Collect the attribute vectors into subsets corresponding to labels.
    Xs <- list()
    for (i in labels) Xs[[i]] <- X[y == i,]
    
    # Now pair up the attribute vector sets and label the first of the
    # pair +1 and the second of the pair -1. Randomly permute the
    # result if necessary.
    k <- 1
    result = list()
    for (p in pairClasses(Xs, labels)) {
        Xpos <- p[[1]]
        Xneg <- p[[2]]
        m1 <- nrow(Xpos)
        m2 <- nrow(Xneg)
        thisX <- rbind(Xpos,Xneg)
        thisY <- c(rep(1, m1), rep(-1, m2))
        m <- m1+m2
        if (randomize) perm <- sample(1:m)
        		  else perm <- (1:m)
        result[[k]] <- list(X = thisX[perm,], y = thisY[perm], 
        					classes = list(p[[3]], p[[4]]))
        k <- k+1
    }
    return(result)
}

#-------------------------------------------------------------------------------
# Check the output of lpSVMLearnMulti to see if any optimizations failed.
# If all is well the result should be 0. Otherwise you get the number of
# failed optimizations.
#-------------------------------------------------------------------------------
checkForFailures <- function(result) {
	count <- 0
	for (r in result) 
		if (r$parameters$status != 0)
			count <- count + 1
	return(count)
}
