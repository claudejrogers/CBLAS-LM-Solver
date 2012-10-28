# CBLAS-LM-Solver

An implementation of the Levenberg-Marquardt algorithm using the cblas library included in the Accelerate framework (Mac OS X).

## Currently supported functions

boltzmann - Boltzmann sigmoid
expdecay  - Exponential decay
gaussian  - Gaussian function
hill      - Hill equation
ic50      - Dose response
mm        - Michaelis-Menten
modsin    - Sine wave


## Adding new functions

This program can be extended to accommodate more functions than those currently included.
Say you want to add the following function:

![model](/images/model.png)

Where

![params](/images/params.png)

The goal of the program is to find values for **x** that minimize the sum of the square residuals between the data and the model.

* Provide a name for the model, say `expfunc`, and add the model name to the `enum Model` in `CblasLMSolver.h`

```c
enum Model {
    boltzmann,
    expdecay,
    expfunc, // newly added
    gaussian,
    hill,
    ic50,
    mm,
    modsin
};
```

* Edit the `switch` statement in the `initWithFileAtPath` function in `CblasLMSolver.m` (there are 4 variable to solve in this model).

```c
switch (model) {
    case boltzmann:
    case expfunc: // newly added
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
        break;
    default:
        NSLog(@"varlen fail.");
        exit(1);
        break;
}
```

* Next, edit the `switch` statement in the `equation` method in this file:

```c
...
    case expfunc:
        m[i] = var[2] * exp(var[0] * xi) + var[3] * exp(var[2] * t);
        break
...
```
Note that the **x** vector of the equation above is represented as the array `var`, and the independent variable _t_ is represented by `xi`.

* Compute the partial derivatives of M with respect to each variable in `var`.

![derivatives](/images/derivatives.png)

* Edit the `switch` statement in `derivatives` method in `CblasLMSolver.m`:

```c
...
    case expfunc:
        d1[i] = var[2] * xi * exp(var[0] * xi);
        d2[i] = var[3] * xi * exp(var[1] * xi);
        d3[i] = exp(var[0] * xi);
        d4[i] = exp(var[1] * xi);
        break;
...
```

* Edit the `switch` statement in `get_jac` method in `CblasLMSolver.m`:

```c
...
    case boltzmann:
    case expfunc: // newly added
    case gaussian:
        jac[i] = -d1[i];
        jac[j] = -d2[i];
        jac[k] = -d3[i];
        jac[l] = -d4[i];
        break;
...
```

* Edit `main` function in `main.m` to include the new function:

```c
if ([modelName isEqualToString:@"boltzmann"]) {
    myModel = boltzmann;
} else if ([modelName isEqualToString:@"expdecay"]) {
    myModel = expdecay;
} else if ([modelName isEqualToString:@"expfunc"]) {     // newly added
    myModel = expfunc;                                   // newly added
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
```

* Edit `usage` function in `main.m`.
