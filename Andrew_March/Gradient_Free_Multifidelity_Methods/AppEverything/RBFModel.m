% RBFModel: Builds a fully linear radial basis function model of the
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
%       krig - updated Kriging model
%       linear - boolean if the model is fully linear, should always be
%       true
%
function [krig,linear]=RBFModel(fn,fid,xk,fxk,krig,theta1,theta2,theta3,theta4,delta,deltamax,pmax,optTheta,fn_params)
deltamax=delta;
D=krig.points;
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
        fhigh=fn(xk+Z(:,i),fid(2),fn_params{:});
        krig.fvals=[krig.fvals;fhigh-fn(xk+Z(:,i),fid(1),fn_params{:})];
        index=[index,length(krig.fvals)];
    end
    linear=true;
end
% Add the new information to the database:
krig.points=D;
% Run AddPoints [use equations 4.3 and 4.4 in Ref 19] to generate the
%       fully linear model:
krig=AddPointsUnc(xk,fxk(2)-fxk(1),D,Y2,krig,index,theta2,theta4*deltamax,pmax,optTheta);
