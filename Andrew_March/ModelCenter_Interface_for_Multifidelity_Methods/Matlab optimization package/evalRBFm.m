% Function evalRBFm
%   Computes the maximum likelihood estimate for the high-fidelity function
%   using the two low-fidelity models.
%
%   Inputs:
%       x - points to evaluate RBF at, points are columnwise
%       fn - handle to the function to use for the low-fidelity estimates
%       fid - fidelity levels to use, fid(1) and fid(2)
%       krigM, krigL - kriging error models
%       fn_param - additional parameters to pass to the low-fidelity
%           functions
%
%   Outputs:
%       fest - maximum likelihood estimate for fhigh(x)
%       var - mean square error estimate for fest(x)
%
function [fest,var]=evalRBFm(x,fn,fid,krigM,krigL,fn_param)
[n,npts]=size(x);
d=2;

% Compute the medium-fidelity information:
f2=zeros(npts,1);
fest=f2;
var2=zeros(npts,1);
var=var2;
Y=krigM.Y2;
ninterps=length(Y(1,:));
for i=1:npts
    point=x(:,i);
    Phi=zeros(1,ninterps);
    for j=1:ninterps
        Phi(j)=phi(norm(point-Y(:,j)),krigM.corr,krigM.params);
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
    f2(i)=fn(point,fid(2),fn_param{:})+[Phi,Pi_t]*krigM.Coeffs;
    [blah1,blah2]=size(krigM.L);
    if(blah1==n+1)
        var2(i)=1;%krigM.sigma_2;
    else
        v=krigM.L\(Phi'); 
        var2(i)=krigM.sigma_2*(phi(norm(point-point),krigM.corr,krigM.params)-v'*v);
    end
end

% Compute the low-fidelity information:
f1=zeros(npts,1);
var1=zeros(npts,1);
Y=krigL.Y2;
ninterps=length(Y(1,:));
for i=1:npts
    point=x(:,i);
    Phi=zeros(1,ninterps);
    for j=1:ninterps
        Phi(j)=phi(norm(point-Y(:,j)),krigL.corr,krigL.params);
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
    f1(i)=fn(point,fid(1),fn_param{:})+[Phi,Pi_t]*krigL.Coeffs;
    if(blah1==n+1)
        var1(i)=1;%krigL.sigma_2;
    else
        v=krigL.L\(Phi'); % check transpose...
        var1(i)=krigL.sigma_2*(phi(norm(point-point),krigL.corr,krigL.params)-v'*v);
    end
    
    % Combine the estimate using the maximum likelihood estimator:
    if(abs(var1(i))+abs(var2(i)) <= 1e-10)
        var(i)=0;
        fest(i)=f1(i);
    else
        fest(i)=f1(i)*(var2(i)/(var2(i)+var1(i)))+f2(i)*(var1(i)/(var2(i)+var1(i)));
        var(i)=1/(1/var1(i)+1/var2(i));
    end
end