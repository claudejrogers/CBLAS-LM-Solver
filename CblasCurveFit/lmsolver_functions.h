//
//  lmsolver_functions.h
//  CblasCurveFit
//
//  Created by Claude Rogers on 7/9/11.
//  Copyright 2011 California Institute of Technology. All rights reserved.
//


/*
 * Find h that solves linear equations of the form Ah = g
 * A is a matrix represented as a row major array of lenght N
 * g and h are vectors
 */
void linalg_solve (double *A, int N, double *g, double *h);

/*
 * Calculates rho
 * fvect and newfvect are vectors of length data
 * h and g are vectors of lenght vars
 * mu is a scalar
 */
double get_rho (double *fvect, double *newfvect, int vars, int data, 
                double mu, double *h, double *g);

/*
 * Calculates mu, the max value of the diagonal element of a
 * a is a matrix represented as a row major array of length N
 */
double get_mu (double *a, int N);

/*
 * Finds values for h that solves (A + mu*I)h = -g
 * a is a matrix represented as a row major array of length N
 * h and g are vectors
 * mu is a scalar
 */
void solve_for_h (double *h, double *a, double *g, double mu, int N);