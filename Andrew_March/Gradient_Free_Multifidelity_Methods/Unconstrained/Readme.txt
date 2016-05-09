Readme.txt

This file discusses all of the files included in this distribution. The 
most important file is AIAA-2010-2912-325.pdf which provides a discussion
of the methods, the underlying algorithm, and references.

%   ********** MAIN OPTIMIZATION ROUTINES ******************

Airfoil_Optimization.m: The main file to optimize an airfoil with two-
    fidelity levels.

Airfoil_Multi.m: The main file to optimize an airfoil with three-fidelity
    levels, 1 high-fidelity, and 2 lower-fidelity.



%   ********** Airfoil Analysis Routines ******************

analyze.m: This is the function that performs the airfoil analysis. It
    selects the appropriate fidelity level to use and then calls one of the
    following functions to actually perform the analysis:

CartMesh.m: This function creates a surface mesh for an airfoil and then
    creates Cart3D input files, and uses the scripts in cart and callCart 
    to run Cart3D. -- These scripts will need to be modified and are
    platform dependent.

shockExp.m: This function estimates the lift and drag for a supersonic 
    airfoil using shock-expansion theory.

thickAirfoil.m: This function uses a linearized supersonic flow model to 
    estimate the lift and drag for a supersonic airfoil. It is referred to
    as the panel method in the paper.

thinAirfoil.m: THis function uses a linearized supersonic flow model on the
    camberline of an airfoil. It is referred to as the camberline model in
    the paper and is intended to be a poor model of the airfoil lift and 
    drag estimates.

    Helper files:
        Oblique.m: Used to solve for the oblique shock angle.
        Prandtl.m: Used to solve for the Prandtl-Meyere shock expansion
            angle.

%   ********** Airfoil Geometry Routines ******************
GenFoilSpline.m: Generates an airfoil from surface spline points.

GenFoilPts.m: Generates an airfoil by linearly interpolating surface
    points.

GenFoilFour.m: Generates an airfoil using Fourier coefficients for the 
    camberline and the thickness.

    Helpter files:
        cosspace.m: Cossine spacing function, similar to logspace or 
            linspace.

%   ********** Algorithm Helper Routines ******************

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
