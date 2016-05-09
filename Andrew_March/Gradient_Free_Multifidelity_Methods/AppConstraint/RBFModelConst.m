% RBFModelConst: Builds a fully linear radial basis function model of the
%   high-fidelity constraint from one low-fidelity function.
%
%   Paramters:
%       fn - requires the function handle so it may add points if necessary
%       xk - current iterate
%       fxk - vector of function values at xk, 1 for each 
%           fidelity level
%       D - All locations at which the high-fidelity function has been
%       evaluated.
%       fvals - function values for all points in D
%       constr - function handle the to constraint (not used)
%       theta1-4 -parameters defined in Ref: 17
%        theta2>0, theta4>=theta3>=1, theta1 (0,1/theta3], deltamax>=delta_k>0,
%       delta - trust region size
%       deltamax - maximum trust region size, no longer used
%       pmax - maximum number of calibration points, at least n+1
%       fn_params - cell with additional parameters to use with fn
%       fid - fidelity levels to use with fn
%       optTheta - use maximum likelihood estimator or not (boolean)
%
%   Output:
%       krig - updated Kriging model
%       linear - boolean if the model is fully linear, should always be
%       true
%
function [Y2,Coeffs,linear,D,fvals,params]=RBFModelConst(fn,xk,fxk,D,fvals,constr,theta1,theta2,theta3,theta4,delta,deltamax,pmax,corr,params,fn_param,fid,optTheta)
% theta2>0, theta4>=theta3>=1, theta1 (0,1/theta3], deltamax>=delta_k>0,
% pmax>=n+1
%
deltamax=delta;
% theta4=1;
% theta4=1;

% Begin what should be the function
[Y2,linear,Z,index]=AffPoints(xk,D,theta1,theta3*delta);
if(~linear)
    % need to add n+1-|Y| points
    [Y2,linear,Z,index]=AffPoints(xk,D,theta1*theta3/theta4,theta4*deltamax);
end
if(~linear)
    Z=Z*delta;
    for i=1:length(Z(1,:))
        Y2=[Y2,xk+Z(:,i)];
        D=[D,xk+Z(:,i)];
%         fvals=[fvals;fn(Y2(:,i+1),fid(2),fn_param{:})-fn(Y2(:,i+1),fid(1),fn_param{:})];
        fvals=[fvals;fn(xk+Z(:,i),fid(2),fn_param{:})-fn(xk+Z(:,i),fid(1),fn_param{:})];
        index=[index,length(fvals)];
    end
    linear=true;
end
% Run AddPoints &  use equations 4.3 and 4.4 to generate the model
[Y2,Coeffs,params]=AddPoints(xk,fxk,D,fvals,Y2,index,theta2,theta4*deltamax,pmax,corr,params,optTheta);