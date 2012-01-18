//
//  main.m
//  CblasCurveFit
//
//  Created by Claude Rogers on 7/9/11.
//  Copyright 2011 California Institute of Technology. All rights reserved.
//

#import <Foundation/Foundation.h>
#import "CblasLMSolver.h"

void usage(const char *argv[]);

int main (int argc, const char * argv[])
{
    NSAutoreleasePool * pool = [[NSAutoreleasePool alloc] init];
    NSUserDefaults *args = [NSUserDefaults standardUserDefaults];
    CblasLMSolver *solver;
    // initialize defaults
    enum Model myModel;
    [args setObject:@"" forKey:@"f"];
    [args setDouble:1.0 forKey:@"v1"];
    [args setDouble:1.0 forKey:@"v2"];
    [args setDouble:1.0 forKey:@"v3"];
    [args setDouble:1.0 forKey:@"v4"];
    // get command line args
    NSString *filepath = [args stringForKey:@"f"];
    if ([filepath isEqualToString:@""]) {
        usage(argv);
    }
    NSString *modelName = [args stringForKey:@"m"];
    if ([modelName isEqualToString:@"boltzmann"]) {
        myModel = boltzmann;
    } else if ([modelName isEqualToString:@"expdecay"]) {
        myModel = expdecay;
    } else if ([modelName isEqualToString:@"gaussian"]) {
        myModel = gaussian;
    } else if ([modelName isEqualToString:@"hill"]) {
        myModel = hill;
    } else if ([modelName isEqualToString:@"ic50"]) {
        myModel = ic50;
    } else if ([modelName isEqualToString:@"mm"]) {
        myModel = mm;
    } else if ([modelName isEqualToString:@"modsin"]) {
        myModel = modsin;
    } else {
        usage(argv);
    }
    double v1 = [args doubleForKey:@"v1"];
    double v2 = [args doubleForKey:@"v2"];
    double v3 = [args doubleForKey:@"v3"];
    double v4 = [args doubleForKey:@"v4"];
    NSArray *initialGuess;
    
    initialGuess = [NSArray arrayWithObjects:[NSNumber numberWithDouble:v1], 
                                             [NSNumber numberWithDouble:v2], 
                                             [NSNumber numberWithDouble:v3],
                                             [NSNumber numberWithDouble:v4],
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

void usage(const char *argv[])
{
    printf("\n"
           "Usage: %s -f <FILENAME> -m <MODEL> [-v1 <var1> -v2 <var2> -v3 "
           "<var3> -v4 <var4>]\n"
           "\n"
           "Fit experimental data to a model.\n"
           "\n"
           "FILENAME: Required. Path to input data.\n"
           "    -f FILENAME    File must contain tab delimited x and y data.\n"
           "MODEL:    Required. The model the data is expected to fit.\n"
           "    -m MODEL       Supported models:\n"
           "                   boltzmann - Boltzmann sigmoid\n"
           "                   expdecay  - Exponential decay\n"
           "                   gaussian  - Gaussian function\n"
           "                   hill      - Hill equation\n"
           "                   ic50      - Dose response\n"
           "                   mm        - Michaelis-Menten\n"
           "                   modsin    - Sine wave\n"
           "INITIAL GUESS: Optional. Provide an initial guess to improve "
           "fitting efficiency.\n"
           "Default values are 1.0.\n"
           "                 boltzmann expdecay  gaussian  hill   ic50   mm\n"
           "   -v1 VALUE     min       a         a         delta  delta  Vmax\n"
           "   -v2 VALUE     max       b         b         ic50   ic50   Km\n"
           "   -v3 VALUE     v50       lambda    mu        hill   hill   N/A\n"
           "   -v4 VALUE     slope     N/A       sigma     N/A    N/A    N/A\n"
           "\n"
           "                 modsin\n"
           "                 a\n"
           "                 b\n"
           "                 c\n"
           "                 N/A\n"
           "\nOUTPUT:\n"
           "Updated values for the variables are displayed for each iteration,\n"
           "along with |f(x)|.\n\n", argv[0]);
    exit(0);
}