% Function evalRBF
%   Computes the maximum likelihood estimate for the high-fidelity function
%   using one low-fidelity model.
%
%   Inputs:
%       x - points to evaluate RBF at, points are columnwise
%       fn - handle to the function to use for the low-fidelity estimates
%       Y - Points used in the Kriging model
%       coeffs - The coefficients of the Kriging model
%       corr - The correlation parameter for the Kriging model
%       param - parameters for the Kriging model
%       fn_param - additional parameters to pass to the low-fidelity
%           functions
%       fid - fidelity levels to use, fid(1) and fid(2)
%
%   Outputs:
%       fs - estimate for fhigh(x)
%
function fs=evalRBF(x,fn,Y,coeffs,corr,param,fn_param,fid)
[n,npts]=size(x);
d=2;
fs=zeros(npts,1);
ninterps=length(Y(1,:));
for i=1:npts
    point=x(:,i);
    Phi=zeros(1,ninterps);
    for j=1:ninterps
        Phi(j)=phi(norm(point-Y(:,j)),corr,param);
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
%     for j=1:ninterps
%         fs(i)=fs(i)+coeffs(j)*phi(norm(point-Y(:,j)),corr,param);
%     end
%     for j=ninterps+1:length(coeffs)
%         if(j==ninterps+1)
%             fs(i)=fs(i)+coeffs(j);
%         else
%             fs(i)=fs(i)+coeffs(j)*x(j-ninterps-1,i);
%         end
%     end
%     Phi
%     Pi_t
    fs(i)=fn(point,fid(1),fn_param{:})+[Phi,Pi_t]*coeffs;
end
