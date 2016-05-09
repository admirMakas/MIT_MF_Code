%  phi=Penalty(x,sigma,fn,Y2,Coeffs,corr,params,constr,fn_param,fid)
%   This function does the quadratic penalty function calculation.
%
%   Inputs:
%       x - design being evaluated.
%       sigma - current penalty parameter value
%       fn - function handle to the objective function
%       Y2 - Points used in the current RBF correction
%       Coeffs - Coefficients of the current RBF model
%       corr - correlation used in RBF model
%       params - parameters to the RBF correlation 
%       constr - function handle to the constraints
%       fn_param - parameters to the function in fn
%       fid - fidelity levels
%
%   Outputs:
%       phi - value of the penalty function at xk
%
function phi=Penalty(x,sigma,fn,Y2,Coeffs,corr,params,constr,fn_param,fid)
[c,ceq]=constr(x,fn_param{:});
m=evalRBF(x,fn,Y2,Coeffs,corr,params,fn_param,fid);
c=c(c>0);
phi=m;
if(~isempty(c))
    phi=phi+sigma*c'*c;
end
if(~isempty(ceq))
    phi=phi+sigma*ceq'*ceq;
end