//
//  main.m
//  CblasCurveFit
//
//  Created by Claude Rogers on 7/9/11.
//  Copyright 2011 California Institute of Technology. All rights reserved.
//

#import <Foundation/Foundation.h>
#import "CblasLMSolver.h"

int main (int argc, const char * argv[])
{
    NSAutoreleasePool * pool = [[NSAutoreleasePool alloc] init];

    CblasLMSolver *solver;
    NSArray *initialGuess;
    NSString *filepath = @"/Users/cjrogers/Desktop/curvefit/cprog/examples/ic50_ex.txt";
    enum Model myModel = ic50;
    
    initialGuess = [NSArray arrayWithObjects:[NSNumber numberWithDouble:1.0], 
                                             [NSNumber numberWithDouble:1.0], 
                                             [NSNumber numberWithDouble:1.0], 
                                             nil];
    
    solver = [[CblasLMSolver alloc] initWithFileAtPath:filepath 
                                                 model: myModel
                                              andGuess:initialGuess];
    [solver levenberg_marquardt];
    
    NSLog(@"Success");
    
    [solver release];
    [pool drain];
    return 0;
}

