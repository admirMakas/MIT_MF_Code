Readme.txt

This file discusses all of the files included in this distribution. The 
most important file is AIAA-2010-2912-325.pdf which provides a discussion
of the methods, the underlying algorithm, and references.

%   ********** MAIN OPTIMIZATION ROUTINES ******************

Example_Optimization.m: The main file to optimize a ModelCenter model using
    the unconstrained optimization presented in AIAA-2010-2912.

Example_Multi.m: The main file to optimize a ModelCenter model using the three
    fidelity (two low-fidelity functions and a single high-fidelity function) 
    presented in the unconstrained optimization paper: AIAA-2010-2912.

%   ********** Airfoil Analysis Routines ******************

ModelCenterWrap2.m: The file the calls the ModelCenter API and runs a ModelCenter
    project. This file will NEED TO BE MODIFIED to run other ModelCenter projects.
    In this file the two low-fidelity models are contained within Matlab and the 
    high-fidelity model is in ModelCenter.


ModelCenterWrap.m: Identical to ModelCenterWrap2.m except in this file all models,
    high- and low-fidelity are contained in ModelCenter.


%  ***********  Other Help Routines  **********************
RBFModel.m: Generates the RBF error model to calibrate a single low-
    fidelity function to a high-fidelity function.

RBFModelm.m: Generates the RBF error models for two low-fidelity functions.

evalRBF.m: (two-fidelity case) Evalautes the low-fidelity function and
 	estimates the high-fidelity function. Returns both the estimate of the 
    high-fidelity function and the estimated mean square error.
    If you get an error about the dimension of krig.L, increase theta1 or theta2.

evalRBFm.m: (three-fidelity case) Evalautes the two low-fidelity functions 
    and computes the maximum likelihood estimate for the high-fidelity 
    functions. Returns both the estimate of the high-fidelity function and 
    the estimated mean square error.
    If you get an error about the dimension of krig.L, increase theta1 or theta2.

AffPoints.m: Finds a basis of affinely independent points within the
    vicinity of the current trust region. This assures the model is fully
    linear.

AddPoints.m: Adds additional points to the RBF error model, in addition to
    the affinely independent basis points. These points are only added if
    they will not cause the RBF model coefficients to be unbounded or have
    a large Hessian. This file includes the maximum likelihood estimate for
    the RBF model and that can/should be changed.


mycons.m: The trust region constraint
