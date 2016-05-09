function [Hdx] = HessVectProduct(x,dx)
% HESSVECT
% Hessian-Vector product subroutine
% Copyright (c) 2007 by Tan Bui
% All Rights Reserved. Massachusetts Institute of Technology

h(1,1) = -400*(x(2)-x(1).^2)+ 800*x(1).^2 + 2.0;
h(1,2) = -400*x(1);
h(2,1) = h(1,2);
h(2,2) = 200;

Hdx = h*dx;
