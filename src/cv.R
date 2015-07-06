#-------------------------------------------------------------------------------
# Make folds for stratified K-fold cross-validation. Assumes 
# multiple classes labelled consecutively from 1 to c unless 
# only two classes are detected with labels +1 and -1. 
#
# Deals gracefully with case of only one label appearing, and with 
# missing classes. Explicitly does NOT deal with leave-one-out. 
#
# Copyright (C) Sean Holden 2014-15.
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# data is a list containing data in X and y.
# Apply the same random permutation to each.
#-------------------------------------------------------------------------------
permuteSets <- function(data) {
	thisy <- data$y
	perm <- sample(1:length(thisy))
	list(X = (data$X)[perm,], y = thisy[perm])
} 

#-------------------------------------------------------------------------------
# Change -1 to 2 for binary case.
#-------------------------------------------------------------------------------
convertBinToMany <- function(y) {
	y[y == -1] <- 2
	return(y)
}

#-------------------------------------------------------------------------------
# Change 2 to -1 for binary case.
# Need to reproduce X to make it easy to use with lapply.
#-------------------------------------------------------------------------------
convertManyToBin <- function(data) {
	thisy <- data$y
	negatives <- (thisy == 2)
	thisy[negatives] <- -1
	return(list(X = data$X, y = thisy))
}

#-------------------------------------------------------------------------------
# Make folds from data. data is a list containing X and y.
#
# Randomizing applies a permutation at the start so that each time 
# you call you'll get a different division. Individual folds are 
# always randomized, so if you call twice with randomize=FALSE you'll get 
# the same examples in each fold, but in a different order.
#-------------------------------------------------------------------------------
makeStratifiedCVFolds <- function(data, K, randomize = TRUE, 
								  verbose = FALSE) {
	
	if (randomize) data <- permuteSets(data)
	y <- data$y
	X <- data$X
	n <- length(y)
	
	# You're going to want at least one example in each fold.
	# We're also not dealing with leave-one out as it's much 
	# more sensible to treat that seperately.
	stopifnot((K > 1) && (K < n))
	
	# Number of classes actually appearing in the data, which is not 
    # necessarily the number of classes for the problem. 
    # The usual experimental setup can reduce the size of the sets. 
    classes <- unique(y)
    numClasses <- length(classes)
    
    # Two classes can either mean a binary problem, or a multiclass 
    # problem without all labels represented.
    binary <- ((numClasses == 2) && (1 %in% classes) && (-1 %in% classes))
	if (binary) {
		y <- convertBinToMany(y)
		classes <- c(1,2)
	}
	maxClass <- max(classes)
	
	# How many examples of each class do you have, how many 
	# should be in each fold, and how many are left over?
	# Seperate out the examples and permute at the same time.
	perFold <- rep(0, maxClass)
	leftOver <- rep(0, maxClass)
	classedX <- list()
	for (i in classes) {
		m <- sum(y==i)
		pF <- m %/% K
		if ((pF < 1) && verbose) 
			cat("\nWarning, there will be a fold with no class",i,"\n")
		perFold[i] <- pF
		leftOver[i] <- m %% K
		classedX[[i]] <- X[y == i,]
	}
	
	# Make the initial folds.
	foldSize <- sum(perFold)
	result = list()
	starts <- rep(1, maxClass)
	for (i in 1:K) {
		newX <- NULL
		newy <- NULL
		for (j in classes) {
			pF <- perFold[j]
			if (pF > 0) {
				start <- starts[j]
				end <- start + pF - 1
				newX <- rbind(newX, classedX[[j]][start:end,])
				newy <- c(newy, rep(j, pF))
	        	starts[j] <- end + 1 
			}
		}
		result[[i]] <- list(X=newX,y=newy)
	}

	# Distribute the leftovers
	fold <- 1
	for (i in classes) {
		toAdd <- leftOver[i]
		if (toAdd > 0) {
			p <- starts[i]
			for (j in 1:toAdd) {
				result[[fold]]$X <- rbind(result[[fold]]$X, classedX[[i]][p,])
				result[[fold]]$y <- c(result[[fold]]$y, i)
				p <- p + 1
				fold <- fold + 1
				if (fold > K) fold <- 1
			}
		}
	}	
	
	# Each fold is mostly in order of class so now we definitely 
	# want to randomize.
	result <- lapply(result, permuteSets)
	if (binary) result <- lapply(result, convertManyToBin)

	return(result)
}

#-------------------------------------------------------------------------------
# Given existing folds, construct the training and test 
# set using the kth fold as the test set.
#-------------------------------------------------------------------------------
cvMakeSets <- function(folds, k) {
	K <- length(folds)
	
	thisX <- NULL
	thisy <- NULL
	for (i in 1:K) {
		if (i != k) {
			thisX <- rbind(thisX, folds[[i]]$X)
			thisy <- c(thisy, folds[[i]]$y)
		}
	}
	return(list(X= thisX, y = thisy, testX=folds[[k]]$X, testy=folds[[k]]$y))
}

