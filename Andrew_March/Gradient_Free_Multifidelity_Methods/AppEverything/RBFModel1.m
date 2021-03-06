% RBFModel1: Builds a fully linear radial basis function model of the
%   high-fidelity constraint from one low-fidelity function.
%
%   Paramters:
%       fn - requires the function handle so it may add points if necessary
%       xk - current iterate
%       fxk - Error between functions at xk
%       D - all locations at which the high-fidelity constraint has been
%           evaluated
%       fvals - function value at those locations.
%       modelnum - (place holder) use: 1
%       theta1-4 -parameters defined in Ref: 17
%        theta2>0, theta4>=theta3>=1, theta1 (0,1/theta3], deltamax>=delta_k>0,
%       delta - trust region size
%       deltamax - maximum trust region size, no longer used
%       pmax - maximum number of calibration points, at least n+1
%       corr - RBF correlation function
%       params - parameters for RBF correlation function
%       fn_params - cell with additional parameters to use with fn
%       fid - fidelity levels to use with fn
%       optTheta - use maximum likelihood estimator or not (boolean)
%
%   Output:
%       krig - updated Kriging model
%       D - All locations at which the high-fidelity constraint has been
%           evaluated
%       fvals - value of the high-fidelity constraint at those locations.
%
function [krig,D,fvals]=RBFModel1(fn,xk,fxk,D,fvals,modelnum,theta1,theta2,theta3,theta4,delta,deltamax,pmax,corr,params,fn_param,fid,optTheta)

deltamax=delta;


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
%         fvals=[fvals,fn(Y2(:,i+1),fid(2),fn_param{:})-fn(Y2(:,i+1),fid(1),fn_param{:})];
        fvals=[fvals,fn(xk+Z(:,i),fid(2),fn_param{:})-fn(xk+Z(:,i),fid(1),fn_param{:})];
        index=[index,length(fvals)];
    end
    linear=true;
end
% Run AddPoints &  use equations 4.3 and 4.4 to generate the model
[Y2,Coeffs,params]=AddPoints(xk,fxk(modelnum),D,fvals(modelnum,:)',Y2,index,theta2,theta4*deltamax,pmax,corr,params,optTheta);
% count=0;
% minpts=4;
% while(length(Y2(1,:))<minpts && count<10 && length(D(1,:))>minpts+1)
%     'yes'
%     deltamax=deltamax*10;
%     [Y2,Coeffs]=AddPoints(xk,fxk,D,fvals,Y2,index,theta2,theta4*deltamax,pmax,corr,params);
%     count=count+1;
% end
% RBFModel is the set of points Y2 and coefficients contained in Coeffs
krig.Y2=Y2;
krig.Coeffs=Coeffs;
krig.params=params;
krig.linear=linear;
krig.corr=corr;