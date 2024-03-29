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


#---------------------------------------------------------------------------
#
# Functions for basic demos of lpSVM material.
#
#---------------------------------------------------------------------------

source("easyData.R")
source("lpSVMMultiP-c.R")

#---------------------------------------------------------------------------
# Make a grid for visualizing decision boundaries.
#---------------------------------------------------------------------------
pointGrid <- function(xr, yr, size) {
    x <- seq(from = xr[1], to = xr[2], by = size)
    y <- seq(from = yr[1], to = yr[2], by = size)
    M <- matrix(0, length(x)*length(y), 2)
    k <- 1
    for (i in 1:length(x)) {
        for (j in 1:length(y)) {
            M[k,1] <- x[i]
            M[k,2] <- y[j]
            k <- k+1
        }
    }
    return(M)
}

#----------------------------------------------------------------------------
# Test for four classes with the lpSVM with probabilistic outputs.
#----------------------------------------------------------------------------
demo4P <- function() {
  d <- easyData4(50,50,50,50)
  X <- d$X
  y <- d$y
  plotData4(d)
  
  r <- lpSVMLearnMultiP(d, 1, 1)
  
  if (checkForFailures(r) !=0)
    cat("Failed",checkForFailures(r),"times.")
  else {
    cat("Classify...\n")
    f <- lpSVMMultiPf(r, y, 1)
    cat("Done...\n")
   
    gr <- pointGrid(range(X[,1]), range(X[,2]), 0.1)
    
    ap<-apply(gr, 1, f)
    classes <- rep(0, length(ap))
    for (i in 1:length(ap)) classes[i] <- ap[[i]]$class
    
    ap1<-gr[classes==1,]
    ap2<-gr[classes==2,]
    ap3<-gr[classes==3,]
    ap4<-gr[classes==4,]
    
    points(ap1, col="red",pch=20,cex=0.1)
    points(ap2, col="blue",pch=20,cex=0.1)
    points(ap3, col="green",pch=20,cex=0.1)
    points(ap4, col="purple",pch=20,cex=0.1)
  }
}
