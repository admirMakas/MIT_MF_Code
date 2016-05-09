% Function evalRBF
%   Computes the maximum likelihood estimate for the high-fidelity function
%   using one low-fidelity model.
%
%   Inputs:
%       x - points to evaluate RBF at, points are columnwise
%       fn - handle to the function to use for the low-fidelity estimates
%       fid - fidelity levels to use, fid(1) and fid(2)
%       krig - kriging error model
%       fn_param - additional parameters to pass to the low-fidelity
%           functions
%
%   Outputs:
%       f2 - estimate for fhigh(x)
%       var2 - mean square error estimate for fest(x)
%
function [f2,var2]=evalRBF(x,fn,fid,krig,fn_param)
[n,npts]=size(x);
d=2;

% Compute the high-fidelity information:
f2=zeros(npts,1);
var2=zeros(npts,1);
Y=krig.Y2;
ninterps=length(Y(1,:));
for i=1:npts
    point=x(:,i);
    Phi=zeros(1,ninterps);
    for j=1:ninterps
        Phi(j)=phi(norm(point-Y(:,j)),krig.corr,krig.params);
    end
    p_hat=nchoosek(n+d-1,n);
    Pi_t=zeros(1,p_hat);
    for j=1:p_hat
        if(j==1)
            Pi_t(j)=1;
        else
            Pi_t(j)=x(j-1,i);
        end
    end
    f2(i)=fn(point,fid(1),fn_param{:})+[Phi,Pi_t]*krig.Coeffs;
    [blah1,blah2]=size(krig.L);
    if(blah1<=n+1)
        var2(i)=1;
    else
        v=krig.L\(Phi'); % An error here indicates a higher theta1+theta2 maybe needed
        var2(i)=krig.sigma_2*(phi(norm(point-point),krig.corr,krig.params)-v'*v);
    end
end
