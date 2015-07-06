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

#------------------------------------------------------------------------
# Kernels for SVMs, and a function to compute the Gram matrix.
#
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


