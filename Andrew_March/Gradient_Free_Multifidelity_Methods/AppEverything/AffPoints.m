% Function AffPoints:
%   Generates a set of affinely indepdendent points within (or around) the
%   current trust region. See Ref 17, Algorithm 4.2 for more detail.
%
%   Inputs:
%       xk - Current iterate
%       D - All high-fidelity sample points
%       theta - theta1 in algorithm definition, minimum null-space
%           projection to be affinely independent, theta (0,1]
%       delta - trust region size
%
%   Outputs:
%       Y2 - Affinely independent vectors (points in D)
%       linear - if there are n+1 vectors in Y2 then true (making the basis
%           fully linear)
%       Z-null space of Y2
%       index - indices of basis vectors from D used in Y2
%
function [Y2,linear,Z,index]=AffPoints(xk,D,theta,delta)
n=length(D(:,1));

Y=zeros(n,1);
Y2=xk;
index=-1;
Z=eye(n);
m=length(D(1,:));
p=round(m*rand(1)+.5);
if(p>1)
    vec=[p:1:m,1:1:p-1];
else
    vec=1:1:m;
end
for j=vec
%     D2(:,j)=D(:,j)-xk;
    if(norm(D(:,j)-xk)<=delta)
        proj_z=Z*inv(Z'*Z)*Z'*(D(:,j)-xk);
        if(norm(proj_z)>=delta*theta)
            Y=[Y,D(:,j)-xk];
            % Update Z
            Z=null(Y'); % transpose??
            Y2=[Y2,D(:,j)];
            index=[index,j];
        end
    end
end
if(length(Y(1,:))==n+1)
    linear=true;
else
    linear=false;
end