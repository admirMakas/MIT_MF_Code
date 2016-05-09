% RBFModel1: Builds a fully linear radial basis function model of the
%   high-fidelity function from one low-fidelity function.
%
%   Paramters:
%       fn - requires the function handle so it may add points if necessary
%       fid - fidelity levels to use with fn
%       xk - current iterate
%       fxk - vector of function values at xk, 1 for each 
%           fidelity level
%       krig - kriging model to be updated
%       theta1-4 -parameters defined in Ref: 17
%        theta2>0, theta4>=theta3>=1, theta1 (0,1/theta3], deltamax>=delta_k>0,
%       delta - trust region size
%       deltamax - maximum trust region size, no longer used
%       pmax - maximum number of calibration points, at least n+1
%       optTheta - use maximum likelihood estimator or not (boolean)
%       fn_params - cell with additional parameters to use with fn
%
%   Output:
%       Y2 - The points used in the Kriging model
%       Coeffs - The coefficients of the radial basis coefficients
%       linear - boolean if the model is fully linear, should always be
%           true
%       D - All locations at which the high-fidelity function has been
%           sampled
%       fvals - The function values at all locations at which the function
%           has been sampled
%       params - The parameters of the radial basis coefficient
%
function [Y2,Coeffs,linear,D,fvals,params]=RBFModel1(fn,xk,fxk,D,fvals,constr,theta1,theta2,theta3,theta4,delta,deltamax,pmax,corr,params,fn_param,fid,opt_theta)
%
deltamax=delta;

% Attempt to find a Affinely independent basis with the trust region
[Y2,linear,Z,index]=AffPoints(xk,D,theta1,theta3*delta);

% If that has failed, the model will not be fully linear, so expand the
% search regain to be theta3*delta
if(~linear)
    % need to add n+1-|Y| points
    [Y2,linear,Z,index]=AffPoints(xk,D,theta1*theta3/theta4,theta4*deltamax);
end

% There is still not n+1 affinely independent points, sample the function
% at the needed locations:
if(~linear)
    Z=Z*delta;
    for i=1:length(Z(1,:))
        Y2=[Y2,xk+Z(:,i)];
        D=[D,xk+Z(:,i)];
        fvals=[fvals;fn(xk+Z(:,i),fid(2),fn_param{:})-fn(xk+Z(:,i),fid(1),fn_param{:})];
        index=[index,length(fvals)];
    end
    linear=true;
end
% Run AddPoints &  use equations 4.3 and 4.4 to generate the model
[Y2,Coeffs,params]=AddPoints(xk,fxk,D,fvals,Y2,index,theta2,theta4*deltamax,pmax,corr,params,opt_theta);
