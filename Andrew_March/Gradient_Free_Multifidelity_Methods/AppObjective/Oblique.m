% function resid=Oblique(beta,theta,M,gamma)
%   Used to solve for the Oblique shock angle beta, as a function of
%   turning angle theta, upstream Mach number M, and ratio of specific
%   heats gamma.
%
function resid=Oblique(beta,theta,M,gamma)
resid=tan(beta-theta)/tan(beta)-...
    (2+(gamma-1)*M^2*(sin(beta))^2)/((gamma+1)*M^2*sin(beta)^2);