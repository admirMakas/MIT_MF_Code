% This function is used for the 2-norm trust region constraint.
function [c,ceq]=TR(xks,xk,delta)
c=norm(xks-xk,2)-delta;
ceq=[];