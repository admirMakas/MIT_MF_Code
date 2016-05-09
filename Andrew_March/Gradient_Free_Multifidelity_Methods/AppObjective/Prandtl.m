% function [v]=Prandtl(M,gamma)
%
% Solves for the Prandtl-Meyer Angle as a function of the Mach number and
% ratio of specific heats gamma.
function [v]=Prandtl(M,gamma)
v=sqrt((gamma+1)/(gamma-1))*atan(sqrt((gamma-1)/(gamma+1)*(M^2-1)))-atan(sqrt(M^2-1));