%  phi=Penalty(x,sigma,fn,Y2,Coeffs,corr,params,constr,fn_param,fid)
%   This function computes the quadratic penalty function for the easy 
%   constraints g & h.
%
%   Inputs:
%       x - design being evaluated.
%       fn - function handle to the objective function
%       sigma - current penalty parameter value
%       fn_param - parameters to the function in fn
%       constrNoAppr - function handle to the constraints
%       constr_param - parameters to the function in constrNoAppr
%
%   Outputs:
%       Pvals - value of the penalty function at xk
%
function Pvals=Penalty(x,fn,fn_param,sigma,constrNoAppr,constr_param)
[n,npts]=size(x);

Pvals=zeros(1,npts);

% add inequality and equality constraints that aren't approx
for i=1:npts
    point=x(:,i);
    [g,h]=constrNoAppr(point,constr_param{:});
    if(~isempty(g) && any(g>0))
        g=g(g>0);
    else
        g=0;
    end
    if(isempty(h))
        h=0;
    end

    fval=fn(point,fn_param{:});

    Pvals(i)=fval+sigma*(h'*h+g'*g);
end

