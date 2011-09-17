//
//  lmsolver_functions.c
//  CblasCurveFit
//
//  Created by Claude Rogers on 7/9/11.
//  Copyright 2011 California Institute of Technology. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <vecLib/cblas.h>
#include <vecLib/clapack.h>
#include "lmsolver_functions.h"

void linalg_solve(double *A, int N, double *g, double *h)
{
    int i, j;
    __CLPK_integer n, nrhs, lda, ldb, info;
    double *a;
    __CLPK_integer *ipiv;
    
    ipiv = (__CLPK_integer *) malloc(N * sizeof(__CLPK_integer));
    a = (double *) malloc(N*N*sizeof(double));
    
    n = lda = ldb = N;
    nrhs = 1;
    
    // clapack expects column major arrays
    for (i = 0; i < N; i++) {
        h[i] = g[i];
        for (j = 0; j < N; j++) {
            a[i*N + j] = A[i + j*N];
        }
    }
    
    dgesv_(&n, &nrhs, a, &lda, ipiv, h, &ldb, &info);
    
    if (info) {
        printf("error solving linear equation\n");
        exit(3);
    }
    
    free((void *) ipiv);
    free((void *) a);
}

double get_rho(double *fvect, double *newfvect, int vars, int data, 
               double mu, double *h, double *g)
{
    int i;
    double dF, dF1, dF2, dL;
    double *c;
    
    c = (double *) malloc(vars * sizeof(double));
    
    dF1 = cblas_ddot(data, fvect, 1, fvect, 1);
    dF2 = cblas_ddot(data, newfvect, 1, newfvect, 1);
    
    dF = 0.5 * dF1 - 0.5 * dF2;
    
    for (i = 0; i < vars; i++)
        c[i] = mu * h[i] - g[i];
    
    dL = cblas_ddot(vars, h, 1, c, 1);
    dL *= 0.5;
    
    free((void *) c);
    
    return dF/dL;
}

double get_mu (double *a, int N)
{
    int i;
    int alen = N * N;
    double max = a[0];
    
    for (i = 0; i < alen; i += (N + 1))
        max = (max < a[i]) ? a[i] : max;
    
    return 1e-3 * max;
}

void solve_for_h (double *h, double *a, double *g, double mu, int N)
{
    int i, j, dl;
    dl = N * N;
    double *na, *ng, *muI;
    
    muI = (double *) malloc(dl * sizeof(double));
    na  = (double *) malloc(dl * sizeof(double));
    ng  = (double *) malloc(N * sizeof(double));
    
    // multiply identity matrix by mu
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            if (i == j) {
                muI[i*N + j] = mu;
            } else {
                muI[i*N + j] = 0.0;
            }
        }
    }
    
    for (i = 0; i < dl; i++) 
        na[i] = a[i] + muI[i];
    for (i = 0; i < N; i++) 
        ng[i] = -g[i];
    
    // solve (A + muI)h = -g for h
    linalg_solve(na, N, ng, h);
    
    free((void *) muI);
    free((void *) na);
    free((void *) ng);
}

