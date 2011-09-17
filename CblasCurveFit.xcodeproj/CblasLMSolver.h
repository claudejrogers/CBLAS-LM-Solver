//
//  CblasLMSolver.h
//  CblasCurveFit
//
//  Created by Claude Rogers on 7/9/11.
//  Copyright 2011 California Institute of Technology. All rights reserved.
//

#import <Foundation/Foundation.h>
#include "lmsolver_functions.h"

enum Model {
    ic50,
    mm
};

@interface CblasLMSolver : NSObject {
    enum Model  model;
    int         varlen;
    int         datalen;
    double      *var;
    double      *x;
    double      *y;
    double      *m;
@private
    double      *param;
    double      *d1;
    double      *d2;
    double      *d3;
    double      *d4;
}

- (id)initWithFileAtPath:(NSString *)filePath 
                   model:(enum Model) mdl 
                andGuess: (NSArray *)varInit;
- (void)equation;
- (void)derivatives;
- (void)get_f:(double *)fvect;
- (void)get_jac:(double *)jac;
- (void)levenberg_marquardt;

@end
