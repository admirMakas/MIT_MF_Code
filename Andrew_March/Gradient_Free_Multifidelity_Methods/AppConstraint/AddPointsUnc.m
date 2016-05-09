% Function AddPointsUnc: Computes the Coefficients for the RBF Model -- uses a
%   Greedy technique to find points that may be used. The algorithm
%   requires n+1 affinely independent points in the trust region vicinity,
%   then adds all other high-fidelity samples to the RBF model provided
%   they don't violate conditions for the model to be fully linear.
%
%   Inputs:
%       xk - current iterate
%       fxk - error at xk between high- and low- fidelity functions
%       D - All high-fidelity sample locations
%       Y2 - Affinely independent vectors within TR
%       krig - current kriging model structure
%       index - indices of those affinely independent sample points (in D
%           and krig structure)
%       theta - conditioning requirement theta2 for RBF matrix, see Ref 17
%           pg 3208
%       delta - theta4*delta, max area around trust region to check for
%           points to add
%       pmax - maximum number of points allowed
%       opt_theta - Boolean, use maximum likelihood estimator for
%       hyperparameters (1) or just use value given in params definition
%       (0)
%           For maximum likelihood estimator the range of hyperparameters
%           is defined within this function...it is highlighted in the
%           code.
%
%   Output: krig now contains a fully linear RBF model, note because of the
%   Greedy selection process, two different calls to this function may
%   result to two different models for the same inputs. -- The order the
%   points are checked to build the model is random.
%
%   
%
function krig = AddPointsUnc(xk,fxk,D,Y2,krig,index,theta,delta,pmax,opt_theta)

fvals=krig.fvals;
index(find(index==-1,1))=length(fvals)+1;
fvals=[fvals;fxk];
Y2_o=Y2;
index_o=index;


if(opt_theta)
    % ********   YOU MAY WANT TO CHANGE THIS ************  %
    ntest=15; % number of hyperparameters to test
    thetagauss=linspace(0.1,5.1,ntest); % vector of possible
    %   hyperparameters to test, can be modified.
    
    % ********   END OF YOU MAY WANT TO CHANGE THIS ************  %
    
    psi=zeros(size(thetagauss));
else
    % If not optimizing hyperparameters, use given values:
    ntest=0;
    psi=1;
    thetagauss=krig.params(1);
end

% Build the RBF Model:
for q=1:ntest+1
    L=[];
    Y2=Y2_o;
    index=index_o;
    if(q==ntest+1)
        [blah,i]=min(psi);
        if(~isinf(blah))
            krig.params(1)=thetagauss(i);
        else
            krig.params(1)=max(thetagauss);
        end
    else
        krig.params(1)=thetagauss(q);
    end
    
    
    % Build initial matrices to get Coefficients of fully linear portion
    n=length(D(:,1));
    d=2;
    npts=length(D(1,:));
    [dim1,dim2]=size(Y2);
    Phi=zeros(dim2,dim2);
    p_hat=nchoosek(n+d-1,n);
    Pi=zeros(p_hat,dim2);

    for i=1:dim2
        for j=1:dim2
            r=norm(Y2(:,i)-Y2(:,j));
            Phi(i,j)=phi(r,krig.corr,krig.params); 
            if(i==1)
                Pi(i,j)=1;
            elseif(i<=p_hat)
                Pi(i,j)=Y2(i-1,j);
            end
        end
    end

    % Attempt to add points to the RBF Model, see Ref 19
    i=1;
    added=false;
    while(i<=npts && length(Y2(1,:))<=pmax)
        if(isempty(find(index==i,1)))
            if(norm(D(:,i)-xk)<=delta && norm(D(:,i)-xk)>1e-12)
                phi_y=zeros(length(Phi(:,1)),1);
                for j=1:length(phi_y)
                    r=norm(Y2(:,j)-D(:,i));
                    phi_y(j)=phi(r,krig.corr,krig.params);
                end
                pi_y=zeros(1,p_hat);
                for j=1:p_hat
                    if(j==1)
                        pi_y(j)=1;
                    else
                        pi_y(j)=D(j-1,i);
                    end
                end
                Piyt=[Pi';pi_y];
                Phiy=[Phi,phi_y; phi_y',phi(0,krig.corr,krig.params)];
                Zy=null(Piyt');
                [Ly,check]=chol(Zy'*Phiy*Zy,'lower');
                if(check~=0)
                    tau=0;
                else
                    tau=Ly(end,end);
                end
                if(tau>theta)
                    Z=Zy;
                    L=Ly;
                    Pi=Piyt';
                    Phi=Phiy;
                    index=[index,i];
                    Y2=[Y2,D(:,i)];
                    added=true;
                end
            end
        end
        i=i+1;
    end
    [dim1,dim2]=size(Pi);
    
    % Equations 4.3 & 4.4 in Ref. 19.
    f=fvals(index);
    if(dim2>n+1 && added)
        [Q,R]=qr(Pi');
    %     w=inv(Z'*Phi*Z)*(Z'*f);
        w=(Z'*Phi*Z)\(Z'*f);
        lambda=Z*w;
        v=R\Q'*(f-Phi*lambda);
        Coeffs=[lambda;v];
    %     norm(Coeffs-Coeffs0)
    else
        M=[Phi,Pi';Pi,zeros(dim1,dim1)];
        Rhs=[f;zeros(p_hat,1)];
        Coeffs=inv(M)*Rhs;
%         disp(['thetagauss w/prob=',num2str(params(1)),'.'])
    end
    
    % Error check:
    if(length(Coeffs)~=length(Y2(1,:))+p_hat)
        Y2
        Coeffs
        Phi
        Pi
        if(dim2>n+1)
            lambda
            v
        else
            Rhs
        end
        error('myToolbox:AddPoints:badVector','Coeffs is the wrong length.')

    end
    % Compute likelihood and mean square error estimates:
    if(q~=ntest+1)
        npts=length(Y2(1,:));
        if(added && dim2>n+1 && check==0)
            sigma_2=1/npts*(f-Pi'*Coeffs(npts+1:end))'*inv(Phi)*(f-Pi'*Coeffs(npts+1:end));
            psi(q)=det(Phi)^(1/npts)*sigma_2;
        else
            psi(q)=Inf;
        end
    else
        npts=length(Y2(1,:));
        sigma_2=1/npts*(f-Pi'*Coeffs(npts+1:end))'*inv(Phi)*(f-Pi'*Coeffs(npts+1:end));
    end

end
% Assign results to kriging model
krig.Coeffs=Coeffs;
krig.Y2=Y2;
% krig.L=L;
[L,check]=chol(Phi,'lower');
krig.L=L;
krig.sigma_2=sigma_2;

% if(length(L)~=npts)
%     disp('AddPoints Failed: May need to increase theta1 and/or theta2');
%     failed=true;
%     % by the Cholesky decomposition, likely theta2 needs to be higher.
%     npts
%     check
%     L
%     tau
% end    