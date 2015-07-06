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
# version does not appear to work (b) this implements multiclass 
# classification in the same way as libSVM and (c) it includes the 
# use of probabilistic outputs.
#
# This version is extended to use two kernels - see comments in the
# corresponding files for binary problems for details.
#
#-----------------------------------------------------------------------------

source("lpSVMP2K-c.R")
source("multiShared.R")

#-------------------------------------------------------------------------------
# Convergence test for the iterative versions of binary2multiP.
#-------------------------------------------------------------------------------
iterDone <- function(Qp, p, t) {
  # The stopping condition essentially requires (a) sum of elements
  # of p is 1 and (b) Qp is a vector with each element equal to b.
  k <- length(p)
  constant <- Qp[1] * rep(1,k)
  diff <- Qp - constant
  return((abs(sum(p) - 1) < t) && (max(abs(diff)) < t))
}

#-------------------------------------------------------------------------------
# Use the method of Wu et al, JMLR 2004, "Probability Estimates ... Pairwise 
# Coupling" to take a collection of probabilities for pairwise classification 
# problems over k classes and turn them into a vector of k probabilities.
#
# k is the number of classes. The parameter r is a list, each element having
# three values: rij = PR(i | i or j, x), i and j.
#
# The non-iterative version is used by default as the iterative one often
# fails to converge.
#-------------------------------------------------------------------------------
binary2multiP <- function(rs, labels, iter = FALSE, fast = TRUE,
                          tolerance = 1e-10, max = 100) {
  
  k <- length(labels)
  
  # As the labels can be any positive integers we need to translate to  
  # 1:k
  index <- function(x) { which(labels == x) }
  
  # Extract the rij values. 
  R <- matrix(0, k, k)
  for (r in rs) {
    rij <- r[[1]]
    i <- index(r[[2]])
    j <- index(r[[3]])
    R[i,j] <- rij
    R[j,i] <- (1 - rij)
  }
  
  # Set up the Q matrix
  Q <- matrix(0, k, k)
  for (i in 1:k) {
    for (j in 1:k) {
      if (i == j) {
        element = 0
        for (s in 1:k) {
          if (s != i) {
            rsi <- R[s,i]
            element <- element + (rsi^2)
          }
        }
        Q[i,j] <- element
      }
      else {
        Q[i,j] <- -R[j,i] * R[i,j] 
      }
    }
  }

  # Usually you're going to do this using the iteration but it's nice
  # have the option of doing a direct solution.
  if (iter) {

    # The original paper has a simple O(k^2) update and a more subtle O(k)
    # version. However the latter seems tricky so at the moment we'll also allow
    # the former. In fact it's tricky because equation (48) appears incorrect.
    # The term with the sum should be 2 delta (Qp)_t

    # Initialize
    p <- rep(1/k, k)
    t <- 1
    total <- 1

    if (fast) {
      
      # Initialize 
      Qp <- Q %*% p
      pQp <- t(p) %*% Qp
  
      # Iterate
      repeat {
        Qtt <- Q[t,t]
        delta <- ((1/Qtt) * (-(Qp[t]) + pQp))
        pbar <- p
        pbar[t] <- p[t] + delta
        pn <- pbar/(1 + delta)
        Qpn <- (Qp + (Q[,t] * delta))
        Qpn <- Qpn * rep((1/(1+delta)), k)
        pnQpn <- (pQp + (2 * delta * Qp[t]) + (Qtt * (delta^2)))/((1 + delta)^2)

        p <- pn
        Qp <- Qpn
        pQp <- pnQpn
        
        if (iterDone(Qp, p, tolerance)) break

        total <- total + 1
        if (total > max) {
          cat("\nbinary2multi: maximum iterations reached.\n")
          break
        }
        t <- t + 1
        if (t > k) {
          t <- 1
          p <- p/(sum(p))    # Not mentioned in paper but seems a good idea.
          Qp <- Q %*% p
          pQp <- t(p) %*% Qp
        }
      }
    }
    # This is the simpler update
    else { 

      # Iterate
      repeat {
        Qtt <- Q[t,t]
        pQp <- t(p) %*% Q %*% p
        sum <- 0
        for (j in 1:k) {
          if (j != t) 
            sum <- sum + (Q[t,j] * p[j])
        }
        p[t] <- (1/Qtt) * (-sum + pQp)
        p <- p/sum(p)
        
        Qp <- Q %*% p
        if (iterDone(Qp, p, tolerance)) break

        total <- total + 1
        if (total > max) {
          cat("\nbinary2multi: maximum iterations reached.\n")
          break
        }
        t <- t + 1
        if (t > k) t <- 1
      }
    }
  }

  # Solve by finding a direct solution.  
  else { 
    S <- cbind(rbind(Q, rep(1,k)), c(rep(1,k),0))
    solution <- solve(S,c(rep(0, k), 1))
    p <- solution[1:k]
  }

  return(p)
}

#-------------------------------------------------------------------------------
# As lpSVMLearnMulti, but using probabilistic outputs.
#
# Here the output is a list with one entry for each binary classifier. The
# format of each entry is:
#
# alpha b status (from training)
# coeffs --- A B status max) (from binary probability estimate)
# data   --- Xbin ybin classes([[1]] label for +1 [[2]] label for -1))
#
# For a description of the parameters see lpSVMlearn in the file
# lpSVM2K-c.R
#-------------------------------------------------------------------------------
lpSVMLearnMultiP <- function(data, X2, C, gamma, gamma2, balanced = FALSE) {
  y <- data$y
  stopifnot(length(unique(y)) > 1)
  
  # Need to join the X and X2 data together to split into binary problems.	
  X <- data$X
  c1 <- ncol(X)
  c2 <- ncol(X2)
  binary <- multiToBinary(list(X = cbind(X, X2), y = y))
  numBinary <- length(binary)
		
  svmLearn <- function(d, x2) { lpSVMLearnP(d, C, gamma, balanced, x2, gamma2) }
  
  result <- list()
  for (i in 1:numBinary) {
  	# Need to split the X data again.
	binAllX <- binary[[i]]$X
	binX2 <- binAllX[,(c1+1):(c1+c2)]
	binary[[i]]$X <- binAllX[,1:c1]  
	binary[[i]]$X2 <- binX2
	currentParameters <- svmLearn(binary[[i]], binX2)
    # The result needs to maintain the training set used as this is needed for 
    # classification.
    result[[i]] <- list(parameters=currentParameters, data=binary[[i]])
    if (currentParameters$status != 0) cat("\nDoh! Non-zero status...\n")
  }
  return(result)
}

#-------------------------------------------------------------------------------
# As lpSVMMulti, but using probabilistic outputs.
#-------------------------------------------------------------------------------
lpSVMMultiP <- function(pList, y, gamma, gamma2, x, x2) {
  numBinary <- length(pList)
  classes <- unique(y)
  numClasses <- length(classes)
  maxClass <- max(classes)
		
  # First make a vector of the probability output by each
  # binary classifier. (The parametersList contains for each
  # binary problem a pair denoting which two classes it addresses.)
  binProbs <- list()	
  for (i in 1:numBinary) {
    p <- pList[[i]]
    prob <- list()
    prob[[1]] <- lpSVMP(p$parameters, p$data, gamma, x, p$data$X2, gamma2, x2)
    prob[[2]] <- p$data$classes[[1]]
    prob[[3]] <- p$data$classes[[2]]
    binProbs[[i]] <- prob
  }

  # Now use the probability values for binary classifiers to get
  # a probability for each class.
  classProbs <- binary2multiP(binProbs, classes)

  # Now see which class has won.
  winners <- classes[classProbs == max(classProbs)]
  
  # Be careful to make sure that the class labels match the 
  # probabilities in the returned distribution.
  dist <- list()
  for (i in 1:numClasses) {
  	dist[[i]] <- list(label=classes[i], probability=classProbs[i])
  }
  
  # Taking just the first winner is dumb but it's also what libSVM does.
  return(list(class=winners[1], distribution=dist))
  
}

#-------------------------------------------------------------------------------
# As lpSVMMultif, but using probabilistic outputs.
#-------------------------------------------------------------------------------
lpSVMMultiPf <- function(plist, d1, y, gamma, gamma2) {
	function(input) { 
		x <- input[1:d1]
		x2 <- input[(d1+1):length(input)]
		lpSVMMultiP(plist, y, gamma, gamma2, x, x2) 
	}
} 
