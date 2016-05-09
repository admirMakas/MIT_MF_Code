%  function [c,ceq]=mycons(xks,xk,delta,constr,fn_param)
%
%  This function takes all of the user-defined constraints g(x) and h(x)
%  and adds the trust region constraint.
%
%  Inputs:
%       xks - xk+sk
%       xk - design iterate
%       delta - trust region size
%       constr - function handle for constraints
%       fn_param - Additional parameters for constraints
%
%   Outputs: (see fmincon for details)
%       c - vector of inequality constraints
%       ceq - vector of equality constraints
%
function [c,ceq]=mycons(xks,xk,delta,constr,fn_param)
[c,ceq]=constr(xks,fn_param{:});
c=[c;norm(xks-xk,2)-delta];