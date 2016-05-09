% function mycons(s,delta)
%   This is the trust region constraint using the 2-norm. Can easily be
%   changed in this file. Of the LB,UB components of fmincon can be
%   modified for the infinity norm (probably more efficient).
%
function [c,ceq]=mycons(s,delta)
c=norm(s,2)-delta;
ceq=[];