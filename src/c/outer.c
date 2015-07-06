/*----------------------------------------------------------------------*/
/* C code for speeding up R computations involving computations with    */
/* kernels. In particular computing the Gram matrix.                    */
/*                                                                      */
/* Linear Programming SVM Transfer Learning                             */
/* https://github.com/ComputationalProteomicsUnit/lpsvm-tl-code         */
/* Copyright (C) 2015  Sean Holden                                      */
/*                                                                      */
/* This program is free software; you can redistribute it and/or modify */
/* it under the terms of the GNU General Public License as published by */
/* the Free Software Foundation; either version 2 of the License, or    */
/* (at your option) any later version.                                  */
/*                                                                      */
/* This program is distributed in the hope that it will be useful,      */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of       */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        */
/* GNU General Public License for more details.                         */
/*                                                                      */
/* You should have received a copy of the GNU General Public License    */
/* along with this program; if not, write to the Free Software          */
/* Foundation, Inc.,                                                    */
/* 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.          */
/*                                                                      */
/*----------------------------------------------------------------------*/

#include<stdio.h>
#include<math.h>
#include<R.h>

/*--------------------------------------------------------------------*/
/* Compute the inner product of two vectors.                          */
/*--------------------------------------------------------------------*/
double inner(double* x1, double* x2, int length) {
  double sum = 0;
  double diff;
  int i;
  for (i = 0; i < length; i++) {
    diff = x1[i] - x2[i];
    sum += diff * diff;
  }
  return sum;
}

/*--------------------------------------------------------------------*/
/* The matrices x and y should have one example per *column* because  */
/* Rs method of representing matrices as vectors then makes the inner */
/* products easier and faster to handle.                              */
/*--------------------------------------------------------------------*/
void makeOuterRBF(double* x, double* y, int* m, int* d, double* gamma, 
		  double* result) {
  int i, j;
  int mval = *m;
  int dval = *d;
  double gval = *gamma;
  for (i = 0; i < mval; i++) 
    for (j = 0; j < mval; j++) 
	result[(j * mval) + i] = 
	  exp(-gval * inner(x + (i * dval), y + (j * dval), dval));
}

/*--------------------------------------------------------------------*/
/* For the Gram matrix we know we have 1s on the diagonal and it's    */
/* symmetric.                                                         */
/*--------------------------------------------------------------------*/
void makeGramRBF(double* x, int* m, int* d, double* gamma,
		 double* result) {
  int i, j;
  double value;
  int mval = *m;
  int dval = *d;
  double gval = *gamma;
  for (i = 0; i < mval; i++) {
    for (j = i; j < mval; j++) {
      if (i == j)
	result[(j * mval) + i] = 1;
      else {
	value = exp(-gval * inner(x + (i * dval), x + (j * dval), dval));
	result[(j * mval) + i] = value;
	result[(i * mval) + j] = value;
      }
    }
  }
}

/*--------------------------------------------------------------------*/
/* Make vector of all Kernel values for a matrix of examples paired   */
/* with a single example.                                             */
/*--------------------------------------------------------------------*/
void makeVectorRBF(double* x, double* newx, int* m, int* d, double* gamma,
		   double* result) { 
  int i;
  int mval = *m;
  int dval = *d;
  double gval = *gamma;
  for(i=0; i < mval; i++) 
    result[i] = exp(-gval * inner(x + (i * dval), newx, dval));
}



