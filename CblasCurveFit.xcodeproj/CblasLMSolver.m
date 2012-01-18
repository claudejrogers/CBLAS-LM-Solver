//
//  CblasLMSolver.m
//  CblasCurveFit
//
//  Created by Claude Rogers on 7/9/11.
//  Copyright 2011 California Institute of Technology. All rights reserved.
//

#define NORM_G 1e-10
#define NORM_H 1e-15
#define MAX_ITER 500
#define LMS_MAX(a, b) ((a) > (b) ? (a):(b))

#include <vecLib/cblas.h>
#import "CblasLMSolver.h"

@implementation CblasLMSolver

- (id) initWithFileAtPath:(NSString *)filePath 
                    model:(enum Model) mdl 
                 andGuess:(NSArray *)varInit
{
    self = [super init];
    if (self) {
        int                 i;
        NSString            *fileContents;
        NSArray             *lines;
        NSArray             *columns;
        NSMutableArray      *xvals;
        NSMutableArray      *yvals;
        NSNumber            *tempX;
        NSNumber            *tempY;
        NSNumberFormatter   *strToNum;
        
        if (![[NSFileManager defaultManager] fileExistsAtPath:filePath]) {
            NSLog(@"File not found at path: %@", filePath);
            exit(0);
        }
        fileContents = [NSString stringWithContentsOfFile:filePath 
                                                 encoding:NSUTF8StringEncoding 
                                                    error:nil];
        lines = [fileContents componentsSeparatedByString:@"\n"];
        strToNum = [[NSNumberFormatter alloc] init];
        xvals = [[NSMutableArray alloc] init];
        yvals = [[NSMutableArray alloc] init];
        for (NSString *line in lines) {
            columns = [line componentsSeparatedByString:@"\t"];
            tempX = [strToNum numberFromString:[columns objectAtIndex:0]];
            tempY = [strToNum numberFromString:[columns objectAtIndex:1]];
            [xvals addObject:tempX];
            [yvals addObject:tempY];
        }
        datalen = (int)[xvals count];
        varlen = 3;
        model = mdl;
        switch (model) {
            case boltzmann:
            case gaussian:
                varlen = 4;
                break;
            case expdecay:
            case hill:
            case ic50:
            case modsin:
                varlen = 3;
                break;
            case mm:
                varlen = 2;
            default:
                NSLog(@"varlen fail.");
                exit(1);
                break;
        }
        
        var   = (double *) malloc(sizeof(double) * varlen);
        param = (double *) malloc(sizeof(double) * varlen);
        
        x  = (double *) malloc(sizeof(double) * datalen);
        y  = (double *) malloc(sizeof(double) * datalen);
        m  = (double *) malloc(sizeof(double) * datalen);
        d1 = (double *) malloc(sizeof(double) * datalen);
        d2 = (double *) malloc(sizeof(double) * datalen);
        d3 = (double *) malloc(sizeof(double) * datalen);
        d4 = (double *) malloc(sizeof(double) * datalen);
                
        for (i = 0; i < varlen; i++) {
            var[i] = [[varInit objectAtIndex:i] doubleValue];
            param[i] = var[i];
        }
        
        for (i = 0; i < datalen; i++) {
            x[i] = [[xvals objectAtIndex:i] doubleValue];
            y[i] = [[yvals objectAtIndex:i] doubleValue];
        }
        
        [strToNum release];
    } 
    return self;
}

- (void)equation
{
    int i;
    for (i = 0; i < datalen; i++) {
        double xi = x[i];
        switch (model) {
            case boltzmann:
                m[i] = var[0] + ((var[1] - var[0]) /
                                 (1 + exp((var[2] - xi)/var[3])));
                break;
            case expdecay:
                m[i] = var[0] + var[1] * exp(-var[2] * xi);
                break;
            case gaussian:
                m[i] = var[0] + var[1] * exp(-(xi - var[2]) * (xi - var[2]) /
                                             (var[3] * var[3]));
                break;
            case hill:
                m[i] = (var[0] / (1 + pow((var[1]/xi), var[2])));
                break;
            case ic50:
                m[i] = (1 - (var[0]/(1 + pow((var[1]/xi), var[2]))));
                break;
            case mm:
                m[i] = ((var[0] * xi) / (var[1] + xi));
                break;
            case modsin:
                m[i] = var[0] * sin(M_PI * (xi - var[1]) / var[2]);
                break;
            default:
                exit(1);
                break;
        }
    }
}

- (void)derivatives
{
    int i;
    for (i = 0; i < datalen; i++) {
        double xi = x[i];
        switch (model) {
            case boltzmann:
                d1[i] = 1 - (1 / (1 + exp((var[2] - xi) / var[3])));
                d2[i] = 1 / (1 + exp((var[2] - xi) / var[3]));
                d3[i] = (((var[0] - var[1]) *
                          exp((var[2] + xi) / var[3])) /
                         (var[3] * pow((exp(var[2] / var[3]) +
                                        exp(xi / var[3])), 2)));
                d4[i] = (((var[1] - var[0]) * (var[2] - xi) *
                         exp((var[2] - xi) / var[3])) / 
                        (var[3] * var[3] *
                         pow((exp((var[2] - xi) / var[3]) + 1), 2)));
                break;
            case expdecay:
                d1[i] = 1.0;
                d2[i] = exp(-var[2]*xi);
                d3[i] = -xi*var[1]*exp(-var[2]*xi);
                break;
            case gaussian:
                d1[i] = 1.0;
                d2[i] = exp(-(xi - var[2])*(xi - var[2])/
                            (var[3] * var[3]));
                d3[i] = (2*((xi - var[2])/(var[3] * var[3]))*var[1]*
                         exp(-(xi - var[2])*(xi - var[2])/(var[3] * var[3])));
                d4[i] = (2*((xi - var[2])*(xi - var[2])/
                            (var[3]*var[3]*var[3]))*var[1]*
                         exp(-(xi - var[2])*(xi - var[2])/(var[3] * var[3])));
                break;
            case hill:
                d1[i] = (1/(1 + pow((var[1]/xi), var[2])));
                d2[i] = ((-var[0] * var[2] * pow(var[1], (var[2] - 1)) *
                          pow(xi, var[2])) / 
                         pow((pow(var[1], var[2]) + pow(xi, var[2])), 2));
                d3[i] = ((var[0]*pow((var[1] * xi), var[2]) * 
                          log(xi/var[1]))/pow((pow(var[1], var[2]) + 
                                               pow(xi, var[2])), 2));
                break;
            case ic50:
                d1[i] = (-(1/(1 + pow((var[1]/xi), var[2]))));
                d2[i] = ((var[0] * var[2] * pow((var[1]/xi), (var[2] - 1))) / 
                         (xi * pow((1 + (pow((var[1]/xi), var[2]))), 2)));
                d3[i] = ((var[0] * pow((var[1]/xi), var[2]) * 
                          log((var[1]/xi))) / 
                         (pow((1 + pow((var[1]/xi), var[2])), 2)));
                break;
            case mm:
                d1[i] = (xi / (var[1] + xi));
                d2[i] = (-(var[0] * xi) / pow((var[1] + xi), 2.0));
                break;
            case modsin:
                d1[i] = sin(M_PI * (xi - var[1]) / var[2]);
                d2[i] = (-var[0] * M_PI * 
                         cos(M_PI * (xi - var[1]) / var[2])) / var[2];
                d3[i] = ((var[0] * M_PI * (var[1] - xi) *
                          cos(M_PI * (xi - var[1]) / var[2])) / pow(var[2], 2));
                break;
            default:
                exit(2);
                break;
        }
    }
}

- (void)get_f:(double *)fvect
{
    [self equation];
    int i;
    for (i = 0; i < datalen; i++) {
        fvect[i] = y[i] - m[i];
    }
}

- (void)get_jac:(double *)jac
{
    [self derivatives];
    int i, j, k, l;
    i = 0;
    j = datalen;
    k = 2 * datalen;
    l = 3 * datalen;
    for (i = 0; i < datalen; i++, j++, k++, l++) {
        switch (model) {
            case boltzmann:
            case gaussian:
                jac[i] = -d1[i];
                jac[j] = -d2[i];
                jac[k] = -d3[i];
                jac[l] = -d4[i];
                break;
            case expdecay:
            case hill:
            case ic50:
            case modsin:
                jac[i] = -d1[i];
                jac[j] = -d2[i];
                jac[k] = -d3[i];
                break;
            case mm:
                jac[i] = -d1[i];
                jac[j] = -d2[i];
                break;
            default:
                exit(3);
                break;
        }
    }
}

- (void)levenberg_marquardt
{
    int i, k, v;
    double mu, rho;
    double *J;
    double *A;
    double *g;
    double *h;
    double *f;
    double *newf;
    
    k = 0;
    v = 2;
    
    g = (double *) malloc(varlen * sizeof(double));
    h = (double *) malloc(varlen * sizeof(double));
    J = (double *) malloc(varlen * datalen * sizeof(double));
    A = (double *) malloc(varlen * varlen * sizeof(double));
    f = (double *) malloc(datalen * sizeof(double));
    newf = (double *) malloc(datalen * sizeof(double));
    
    [self get_f:f];
    [self get_jac:J];
    
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, 
                varlen, varlen, datalen, 1.0, J, 
                datalen, J, datalen, 0.0, A, varlen);
    cblas_dgemv(CblasRowMajor, CblasNoTrans, varlen, datalen, 
                1.0, J, datalen, f, 1, 0.0, g, 1);
    mu = get_mu(A, varlen);
    
    while ((fabs(g[cblas_idamax(varlen, g, 1)]) >= NORM_G) && (k < MAX_ITER)) {
        k++;
        
        solve_for_h(h, A, g, mu, varlen);
        
        for (i = 0; i < varlen; i++) {
            param[i] = var[i];
            var[i] += h[i];
        }
        if ((cblas_dnrm2(varlen, h, 1)) <= NORM_H) {
            printf("var converged\n");
            break;
        }
        
        [self get_f:newf];
        
        rho = get_rho(f, newf, varlen, datalen, mu, h, g);
        
        if (rho > 0) {
            for (i = 0; i < datalen; i++) 
                f[i] = newf[i];
            [self get_jac:J];
            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, 
                        varlen, varlen, datalen, 1.0, J, 
                        datalen, J, datalen, 0.0, A, varlen);
            cblas_dgemv(CblasRowMajor, CblasNoTrans, varlen, datalen, 
                        1.0, J, datalen, f, 1, 0.0, g, 1);
            mu *= LMS_MAX((1.0/3.0), (1 - pow((2*rho - 1), 3.0)));
            v = 2;
        } else {
            for (i = 0; i < varlen; i++) 
                var[i] = param[i];
            mu *= v;
            v *= 2;
        }
        printf("iter %2d: var = ", k);
        for (i = 0; i < varlen; i++) {
            printf(" %f", var[i]);
        }
        printf(", |f(x)| = %g\n", cblas_dnrm2(varlen, f, 1));
    }
    
    free((void *) J);
    free((void *) A);
    free((void *) g);
    free((void *) h);
    free((void *) f);
    free((void *) newf);
}

- (void)dealloc
{
    free((void *) var);
    free((void *) param);
    free((void *) x);
    free((void *) y);
    free((void *) m);
    free((void *) d1);
    free((void *) d2);
    free((void *) d3);
    free((void *) d4);
    [super dealloc];
}

@end
