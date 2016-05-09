function [f,g] = TestFunction(x)
% TESTFUNCTIION
% Test function for STIRNCG.m
% Copyright (c) 2007 by Tan Bui
% All Rights Reserved. Massachusetts Institute of Technology

% %% Rosenbrock function
 f = 100*(x(2)-x(1).^2).^2 + (1-x(1)).^2;
 if nargout > 1
     g(1,1) = -400*(x(2)-x(1).^2).*x(1) - 2*(1-x(1));
     g(2,1) = 200*(x(2)-x(1).^2);
 end
