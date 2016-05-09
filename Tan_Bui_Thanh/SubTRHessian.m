function [Hd] = SubTRHessian(zc,d,varargin)
% Tan Bui: Jan 7, 2007
% Hessian-vector product for the Lagrangian
    
load SubTRdata, load SubTRHessdata

H = zeros(3,3);
H(1,1) = a(1)+2*z(3)*c(1); H(1,2) = a(2)+2*z(3)*c(2);
H(1,3) = 2*(z(1)*c(1)+z(2)*c(2));
H(2,2) = a(3)+2*z(3)*c(3); H(2,3) = 2*(z(1)*c(2)+z(2)*c(3));
H(2,1) = H(1,2); H(3,1) = H(1,3); H(3,2) = H(2,3);

dz2 = -H(3,1)*d/H(3,2);
dz3 = -(H(2,1)*d+H(2,2)*dz2)/H(2,3);

Hd = H(1,1)*d+H(1,2)*dz2+H(1,3)*dz3;