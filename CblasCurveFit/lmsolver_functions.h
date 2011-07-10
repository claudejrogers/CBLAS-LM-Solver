//
//  lmsolver_functions.h
//  CblasCurveFit
//
//  Created by Claude Rogers on 7/9/11.
//  Copyright 2011 California Institute of Technology. All rights reserved.
//

void linalg_solve (double *A, int N, double *g, double *h);
double get_rho (double *fvect, double *newfvect, int vars, int data, 
                double mu, double *h, double *g);
double get_mu (double *a, int N);
void solve_for_h (double *h, double *a, double *g, double mu, int N);