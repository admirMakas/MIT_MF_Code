clear all
% $$$ STIRNCG  Copyright and License information for STIRNCG  
% $$$     STIRNCG is the subspace trust region interior reflective Newton-CG method
% $$$
% $$$     The following information is a copy of the License file in the STIRNCG
% $$$     distribution.
% $$$  
% $$$  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $$$  
% $$$   STIRNCG, Copyright (c) 2007 by Tan Bui
% $$$   All Rights Reserved. Massachusetts Institute of Technology
% $$$   
% $$$   STIRNCG License:
% $$$   
% $$$       Your use or distribution of STIRNCG or any modified version of
% $$$       STIRNCG implies that you agree to this License.
% $$$   
% $$$       THIS MATERIAL IS PROVIDED AS IS, WITH ABSOLUTELY NO WARRANTY
% $$$       EXPRESSED OR IMPLIED.  ANY USE IS AT YOUR OWN RISK.
% $$$   
% $$$       Permission is hereby granted to use or copy this program, provided
% $$$       that the Copyright, this License, and the Availability of the original
% $$$       version is retained on all copies.  User documentation of any code that
% $$$       uses STIRNCG or any modified version of STIRNCG code must cite the
% $$$       Copyright, this License, the Availability note, and "Used by permission."
% $$$       Permission to modify the code and to distribute modified code is granted,
% $$$       provided the Copyright, this License, and the Availability note are
% $$$       retained, and a notice that the code was modified is included.  This
% $$$       software was developed with support from the National Science Foundation,
% $$$       and is provided to you free of charge.
% $$$   
% $$$   Availability:
% $$$   
% $$$       web.mit.edu/tanbui/Public/TRNCGcode/
% $$$  
% $$$   Main reference:
% $$$
% $$$   1. A Subspace, Interior, and Conjugate Gradient Method for Large-Scale 
% $$$   Bound-Constrained Minimization Problems. SIAM Journal on Scientific
% $$$   Computing, Volume 21 ,  Issue 1  (Aug.-Sept. 1999), Pages: 1 - 23
% $$$ 
% $$$   2. Numerical Optimization, by Nocedal and Wright
% $$$
% $$$   Bug report:
% $$$   tanbui@alum.mit.edu or tansweet@gmail.com
% $$$
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%--------Options for the optimization solvers-----------------------------
option.OptMethod = 'STIRNCG';   % Optimization solver name
option.Scaled = 1;              % (1 or 0) Scale when meet the
                                % bound constraints  or not

option.TolF  = 1.e-12;          % tolerance for function value
option.TolX  = 1.e-10;          % tolerance for optimization step
option.TolGrad = 1.e-8;         % tolerance for the scaled gradient

option.MaxIter = 100;           % maximum number of iterations

option.AffineScaling = 'huuIa'; %'CL' Coleman and Li scaling
                                %'CLModified' Modified CL so that
                                % the maximum TR radius is 1
                                % 'Tan' my scaling
                                % 'huuIa': please read my thesis
                                % and its references
                                % the maximum TR radius is 1 and
                                % many more

option.ReflectionFlag = 1;      % 1: Single reflection, 2: reflect all

option.pcv = @No_PC;          % Preconditioning function: Please
                                % put your preconditioner
                                % here. Otherwise the default with
                                % no preconditioning is used

option.HessVect = @HessVectProduct;     % put your Hessian-vector product
                                % subroutine here

%%--------------------------------------------------------------------------

%% File to write for STIRNCG
if strcmpi(option.OptMethod,'STIRNCG')
    option.SubTRFile = fopen('SubTRHistory','w'); 
end


% Initial guess 
x0 = [0.5;0.5];
% lower bound
xlow = [-inf;-inf]; 
% upper bound
xup = [inf;inf];

% The solver
% xopt:     optimal solution
% fopt:     optimal function value
% exitflag: exit flag
[xopt, fopt, exitflag] = TRBoundOptimizationSolver(x0, @TestFunction, xlow, xup,option);

% optimization message output
OptimizationMessage(exitflag);

% close file if STIRNCG is used
if strcmpi(option.OptMethod,'STIRNCG')
  fclose(option.SubTRFile);
end


