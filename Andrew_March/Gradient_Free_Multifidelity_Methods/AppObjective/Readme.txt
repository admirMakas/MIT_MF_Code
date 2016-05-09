Readme.txt

This file discusses all of the files included in this distribution. The 
most important file is MDAO_Constrained.pdf which provides a discussion
of the methods, the underlying algorithm, and references.

%   ********** MAIN OPTIMIZATION ROUTINES ******************

AppObjective.m: The main file to optimize an airfoil (or other function)
    with two-fidelity levels. The high-fidelity objective function is
    approximated, but the constraints are handled without multifidelity
    techniques, Algorithm 1 the included paper.


%   ********** Airfoil Analysis Routines ******************

analyze.m: This is the function that performs the airfoil analysis. It
    selects the appropriate fidelity level to use and then calls one of the
    following functions to actually perform the analysis:


shockExp.m: This function estimates the lift and drag for a supersonic 
    airfoil using shock-expansion theory.

thickAirfoil.m: This function uses a linearized supersonic flow model to 
    estimate the lift and drag for a supersonic airfoil. It is referred to
    as the panel method in the paper.

thinAirfoil.m: This function uses a linearized supersonic flow model on the
    camberline of an airfoil. It is referred to as the camberline model in
    the paper and is intended to be a poor model of the airfoil lift and 
    drag estimates.

myAFconst.m: This is the function that handles the geometric constraints on
    the airfoil. For example, 5% t/c and all thickness positive.

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

RBFModel1.m: Generates the RBF error model to calibrate a single low-
    fidelity function to a high-fidelity function.


evalRBF.m: Evalautes the low-fidelity function and
 	estimates the high-fidelity function. Returns both the estimate of the 
    high-fidelity function.


AffPoints.m: Finds a basis of affinely independent points within the
    vicinity of the current trust region. This assures the model is fully
    linear.

AddPoints.m: Adds additional points to the RBF error model, in addition to
    the affinely independent basis points. These points are only added if
    they will not cause the RBF model coefficients to be unbounded or have
    a large Hessian. This file includes the maximum likelihood estimate for
    the RBF model and that can/should be changed.

Penalty.m: This function computes the quadratic penalty function value.

mycons.m: The function that appends the trust region constraint to all of the
   user-defined constraints.

TR.m: The trust region constraint.