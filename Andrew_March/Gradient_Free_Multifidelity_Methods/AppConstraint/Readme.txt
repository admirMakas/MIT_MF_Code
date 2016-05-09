Readme.txt

OPTIMIZATION OF A SINGLE-FIDELITY OBJECTIVE -- SUBJECT TO A MULTIFIDELITY
CONSTRAINT.

This file discusses all of the files included in this distribution. The 
most important file is MDAO_Constrained.pdf which provides a discussion
of the methods, the underlying algorithm, and references.

%   ********** MAIN OPTIMIZATION ROUTINES ******************

AppConstraint.m: The main file to optimize an airfoil (or other function)
    with two-fidelity levels. The high-fidelity objective function is
    handled directly, but the constraints are approximated with 
    multifidelity techniques, Algorithm 2 the included paper.

ConstSolver.m: This is step 0 in Algorithm 2 in the paper. It performs the
    constrained Algorithm 1 to find a point that is feasible with respect 
    to all of the constraints. The difference between this algorithm and 
    Algorithm 1 in the paper is the stopping condition, any feasible value
    of the high-fidelity constraint that is negative has been added, and
    the parameter "d" described in this constrained paper used to 
    assure the constraint is bounded from below.


unConst.m: This is step 0 in Algorithm 2 in the paper. It performs the
    unconstrained algorithm in AIAA 2010-2912 to find a point that is 
    feasible with respect to only the high-fidelity constraints. The 
    difference between this algorithm and the algorithm in AIAA 2010-2912
    is the stopping condition, any feasible value of the high-fidelity 
    constraint that is negative has been added, and the a parameter "d"
    described in the consstrained paper used to assure the constraint is
    bounded from below.


%  *** Important Auxiliary Functions ******

SurrConst.m: Is the function that computes the value of the surrogate 
    constraint m-bar, all of the unapproximated constraints g and
    h, and if delta is a positive vale the trust region constraint.

Penalty.m: The function the computes the merit function consisting of 
    only the easy constraints g & h, and the objective function.


%   ********** Airfoil Analysis Routines ******************

Dmax.m: This is the function that performs the airfoil analysis. It
    selects the appropriate fidelity level to use and then calls one of the
    following functions below to perform the analysis.

LDest.m: This is the function that performs the airfoil analysis. It
    selects the appropriate fidelity level to use and then calls one of the
    following functions to actually perform the analysis:

shockExp.m: This function estimates the lift and drag for a supersonic 
    airfoil using shock-expansion theory.

thickAirfoil.m: This function uses a linearized supersonic flow model to 
    estimate the lift and drag for a supersonic airfoil. It is referred to
    as the panel method in the paper.

thinAirfoil.m: THis function uses a linearized supersonic flow model on the
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

%%%%%  These Routines are used to find a feasible starting point, if there
      are no easy constraints g, and h: 

  RBFModel.m: Generates the RBF error model to calibrate a single low-
    fidelity function to a high-fidelity function.


  evalRBF.m: Evalautes the low-fidelity constraint and
 	estimates the high-fidelity constraint value. Note, this includes the
	bounded from below modification, and "d" is currently set to 1 in the
	function.
  
%%%%%  These Routines are used to find a feasible starting point, with
	easy constraints g, and h: 

  RBFModelConst.m: Generates the RBF error model to calibrate a single low-
    fidelity function to a high-fidelity function.


  evalRBFConst.m: Evalautes the low-fidelity constraint and
 	estimates the high-fidelity constraint value. Note, this includes the
	bounded from below modification, and "d" is currently set to 1 in the
	function.
  PenaltyConst.m: This function computes the quadratic penalty function value. 
	Note, this formulation includes the parameter d which is set in the function 
	itself, currently 1.

%%%% Other helper functions:
 
AffPoints.m: Finds a basis of affinely independent points within the
    vicinity of the current trust region. This assures the model is fully
    linear.

AddPoints.m: Adds additional points to the RBF error model, in addition to
    the affinely independent basis points. These points are only added if
    they will not cause the RBF model coefficients to be unbounded or have
    a large Hessian. This file includes the maximum likelihood estimate for
    the RBF model and that can/should be changed.

AddPointsUnc.m: The same function as AddPoints.m, but the unconstrained 
   method uses a different data storage method than the constrained. The
   algorithms are otherwise identical.

mycons.m: The trust region constraint.

myconsConst: This function combines the unapproximated constraints and the 
   trust-region constraint.
