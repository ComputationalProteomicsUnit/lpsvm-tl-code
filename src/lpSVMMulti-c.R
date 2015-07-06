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
# Linear programming SVM as described in Wu and Dietterich, "Improving SVM 
# accuracy by Training on Auxilliary Data Sources".
#
# The only difference is that (a) this enforces xi as the unconstrained 
# version does not appear to work and (b) this implements multiclass 
# classification in the same way as libSVM.
#
#-----------------------------------------------------------------------------

source("lpSVM-c.R")
source("multiShared.R")

#-------------------------------------------------------------------------------
# Linear programming SVM learning for multiple-class 
# problems. Uses the binary version lpSVMLearn. We split the
# multiple-class problem into a collection of binary problems
# and train a binary classifier for each.
#
# y should specify classes as positive integers.
#-------------------------------------------------------------------------------
lpSVMLearnMulti <- function(data, C, gamma, balanced = FALSE) {
	stopifnot(length(unique(data$y)) > 1)
	binary <- multiToBinary(data)
	numBinary <- length(binary)
	
	svmLearn <- function(d) { lpSVMLearn(d, C, gamma, balanced) }
	
	result <- list()
	for (i in 1:numBinary) {
		# The result needs to maintain the training set used as this is
		# needed for classification.
		currentParameters <- svmLearn(binary[[i]])
		result[[i]] <- list(parameters = currentParameters, 
							data=binary[[i]])
    }
	return(result)
}

#-------------------------------------------------------------------------------
# Given the parameters from training a multiclass lpSVM with lpSVMLearnMulti
# classify a new example x.
#-------------------------------------------------------------------------------
lpSVMMulti <- function(pList, y, gamma, x) {
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
			out <- lpSVM(p$parameters, p$data, gamma, x) 
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
# As lpSVMMulti but returns a function for classifying new vectors directly.
#-------------------------------------------------------------------------------
lpSVMMultif <- function(plist, y, gamma) {
	function(x) { lpSVMMulti(plist, y, gamma, x) }
} 
