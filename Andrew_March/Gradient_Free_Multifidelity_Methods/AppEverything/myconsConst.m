% function [c,ceq]=myconsConst(xks,xk,delta,constr,constr_param)
%   This function combines the unapproximated constraints and the
%   trust-region constraint.
%
%   % Inputs:
%       xks - design iterate xk+step sk
%       xk - design iterate xk
%       delta - trust region size
%       constr - function handle to unapproximated constraints
%       constr_param - parameters needed for constraint function
%
%   Outputs:
%       c - inequality constraints, g+trust region(fmincon format)
%       ceq - equality constraints, h (fmincon format)
function [c,ceq]=myconsConst(xks,xk,delta,constr,constr_param)

[c,ceq]=constr(xks,constr_param{:});
c=[c;norm(xks-xk,2)-delta];
