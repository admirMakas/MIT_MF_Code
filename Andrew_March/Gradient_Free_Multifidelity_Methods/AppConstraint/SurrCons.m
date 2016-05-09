% [c,ceq]=SurrCons(x,s,delta,constr,otherConstrs,constrNoAppr,krigs,constr_param,fid)
%
%   SurrCons is the constraint function that computes the value of the
%   surrogate constraint m-bar, all of the unapproximated constraints g and
%   h, and if delta is a positive vale the trust region constraint.
%
%   Note: This function function pads the value of m-bar and returns 
%   m-bar+min(c_high(xk),delta/10) and this can and should be modified
%   depending on the optimization routine being used. Also, the value of
%   c_high(xk) is currently read from the workspace, if the optimization
%   routine is made a function this will need to be modified.
%
%   Inputs:
%       x - current design iterate xk
%       s - size of step xk
%       delta - trust region size, or a negative number. If delta>=0 the
%           trust region constraint is added to the bottom of the inequality
%           constraint vector, if delta<0 then it is not added.
%       constr - function handle to the low-fidelity constraint
%       otherConstrs - Boolean, are there other constraints, g & h?
%       constrNoAppr - Function handle to the unapproximated constraints,
%           only used if otherConstrs=true
%       krigs - RBF model for error between high- and low-fidelity
%           constraints
%       constr_param - additional parameters passed to the constraint
%           functions
%
%   Outputs:  
%       c - Inequality constraints, m-bar, g, and trust-region (fmincon
%           format)
%       ceq - Equality constraints, h
%
%
function [c,ceq]=SurrCons(x,s,delta,constr,otherConstrs,constrNoAppr,krigs,constr_param,fid)
nconstrs=length(krigs);
[n,npts]=size(x); % npts must be 1
point=x+s;
d=2;
fs=zeros(nconstrs,npts);
for k=1:nconstrs
    Y=krigs(k).Y2;
    coeffs=krigs(k).Coeffs;
    corr=krigs(k).corr;
    param=krigs(k).params;
    ninterps=length(Y(1,:));
    for i=1:npts
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

        clows=constr(point,fid(1),constr_param{:});
        fs(k,i)=clows(k)+[Phi,Pi_t]*coeffs;
    end
    
end

fs=fs+min(delta/10,evalin('base','c2')); % This can and should be changed 
%           depending on what algorithm is being used in the optimization
%           problem. For an interior point method, just keep fs=fs; If the
%           solver tolerance is poor a larger value may be required. This
%           currently pulls the value of c2 from the workspace...if c2 is
%           not in the workspace either remove it from here or pass it as a
%           parameter.

if(otherConstrs)
    [c0,ceq0]=constrNoAppr(point,constr_param{:});

    if(delta>=0)
        c=[fs;c0;norm(s)-delta];
    else
        c=[fs;c0];
    end
    ceq=ceq0;
else
    if(delta>=0)
        c=[fs; norm(s)-delta];
    else
        c=fs;
    end
    ceq=[];
end