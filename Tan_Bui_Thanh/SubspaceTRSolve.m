function [] = SubspaceTRSolve(a,b,c,posdef,Delta, SubTRFile)
% SUBSPACETRSOLVE
% solving the subspace trust region problem
% Copyright (c) 2007 by Tan Bui
% All Rights Reserved. Massachusetts Institute of Technology
    
if posdef == 1,
  n = length(b);
  
  if n == 2,
    A = [a(1),a(2); a(2), a(3)];
    
% $$$     cn = rcond(A);
% $$$     if (abs(cn) < eps), keyboard, end
    
    alpha = -(A\b);
    temp = alpha(1)^2*c(1)+2*alpha(1)*alpha(2)*c(2)+alpha(2)^2* ...
           c(3);
    
    if (temp > (Delta^2))
      % need to solve a equally-constraint optimization problem
      [alpha] = STRNCG(a,b,c,Delta,SubTRFile);
    else
      z = alpha;
      save SubTRSolution z
    end
  elseif n == 1,
    alpha = -b/a;
    if ((alpha^2*c)>Delta^2),
      z1 = Delta/sqrt(c);  z2 = -z1;
      m1 = z1*b+0.5*z1^2*a; m2 = z2*b+0.5*z2^2*a;
      if (m1 > m2)
        z = z2;
      else
        z = z1;
      end
      save SubTRSolution z
    else
      z = alpha;
      save SubTRSolution z
    end
  else
    n
    error('n must be 1 or 2')
  end
  
elseif posdef == 0,
  n = length(b);
  if (n == 1),
    z1 = Delta/sqrt(c);  z2 = -z1;
    m1 = z1*b+0.5*z1^2*a; m2 = z2*b+0.5*z2^2*a;
    if (m1 > m2)
      z = z2;
    else
      z = z1;
    end
    save SubTRSolution z
  elseif (n == 2)
    alpha = STRNCG(a,b,c,Delta,SubTRFile);
  else
    n
    error('n must be 1 or 2')
  end
else
  posdef
  error('posdef must be 1 or 0')
end

%----------------------------------------------------------

function alpha = STRNCG(a,b,c,Delta,SubTRFile)
% Solve the subspace trust region using the TRNCG method for the
% Lagrangian 
    
option.TolF  = 1.e-15;
option.TolX  = 1.e-10;
option.TolGrad = 1.e-10;
option.MaxIter = 100;
option.HessVect = @SubTRHessian;
% Preconditioner
option.pcv = '';%@HessVectTest;
option.SubTRFile = SubTRFile;

save SubTRdata a b c Delta

% Compute initial guess for alpha
alpha0 = 0;
% $$$     alpha0(2) = sqrt(Delta^2/c(3));
% $$$     alpha0(3) = -1/2;


[alpha,fc,exitflag,histout,XC] = TRNCG(alpha0,@SubTRCost,option);
