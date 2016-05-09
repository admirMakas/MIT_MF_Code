% function [xk,fk,points,fvals]=ConstSolver(fn,points,fvals,corr,params,fid,fn_param,constr,optTheta)
%
%   This function is used to find a feasible point for the multifidelity
%   constraint optimization, essentially Algorithm 1 in theory paper with
%   additional stopping criteria and the parameter "d" added to the merit
%   functions to assure the value of the constraint is bounded from below.
%
%   Inputs:
%       fn - Handle the function being minimized (the constraint)
%       points - All of the locations the high-fidelity constraint has been
%           evaluated.
%       fvals - The high-fidelity constraint value at those locations
%       corr - RBF function correlation
%       params - RBF function parameters
%       fid - fidelity levels at which the constraints can be evaluated
%       fn_param - parameters that will be used calling the constraints
%       constr - function handle to the unapproximated constraints, g & h
%       optTheta - Boolean, use maximum likelihood estimator for RBF
%       calibration?
%
%   Outputs:
%       xk - a feasible design iterate
%       fk - the high-fidelity constraint value at that point
%       points - all locations at which the high-fidelity constraint has
%           been evaluated
%       fvals - the value of the high-fidelity constraint at all of those
%       locations.
%
function [xk,fk,points,fvals]=ConstSolver(fn,points,fvals,corr,params,fid,fn_param,constr,optTheta)

disp('Minimizing High-Fidelity Constraint Violation')

xk=points(:,1);
flow=fn(points(:,1)',fid(1),fn_param{:});
fk=fvals(1)+flow;
mk=fk;
n=length(xk);


% Algorithm Parameters:    ***********************************************
%       See MDAO paper for details
pmax=20;                    % maximum number of calibration points
eta0=0.25;                  % trust region contraction criteria
eta1=.75;                   % trustregion expansion crtiteria
eta2=2;                     % trust region expansion criteria
gamma0=1/2;                 % trust region contraction ratio
gamma1=2;                   % trust region expansion ratio
delta0=max(1,norm(xk,inf)); % Initial trust region size
deltamax=10^3*delta0;       % Maximum trust region size
epsilon=5*10^(-4);          % Convergence tolerance
a=10^-4;                    % Sufficient decrease condition
theta1=10^-3;               % RBF Model Construction Parameter, see Ref 20
theta2=10^-4;               % RBF Model Construction Parameter, see Ref 20
theta3=10;                  % RBF Model Construction Parameter, see Ref 20
theta4=max(sqrt(n),10);     % RBF Model Construction Parameter, see Ref 20
maxiter=200;                % maximum number of iterations
maxsubiter=200;             % maximum number of iterations on TR subproblem 
fd=10^-6;                   % Finite difference length
% sigma                     % Penalty parameter defined in main body,
                            %   growth rate can be modified


iter=0;
delta=delta0;
epsOpt=epsilon/10;

% Build initial surrogate model for the objective function:
[Y2,Coeffs,linear,points,fvals,params]=RBFModelConst(fn,xk,fk-flow,points,fvals,constr,theta1,theta2,theta3,theta4,delta,deltamax,pmax,corr,params,fn_param,fid,optTheta);


% Compute the initial Lagrange multipliers:
dm_dx=zeros(n,1);
mk=evalRBFConst(xk,fn,Y2,Coeffs,corr,params,fn_param,fid);
[c,ceq]=constr(xk,fn_param{:});
ncon=length(c);
nconE=length(ceq);
dc_dx=zeros(n,ncon);
dceq_dx=zeros(n,nconE);
ConMat=zeros(n,ncon+nconE);
for i=1:n
    dx=zeros(n,1);
    dx(i)=fd;
    dm_dx(i)=(evalRBFConst(xk+dx,fn,Y2,Coeffs,corr,params,fn_param,fid)-mk)/(fd);
    [c2,ceq2]=constr(xk+dx,fn_param{:});
    if(ncon>=1)
        dc_dx(i,:)=(c2'-c')/fd;
    end
    if(nconE>=1)
        dceq_dx(i,:)=(ceq2'-ceq')/fd;
    end
end

indices=[];
c2=0;
for i=1:ncon
    if(c(i)>=-epsOpt && dc_dx(:,i)'*dm_dx<=0)
        indices=[indices,i];
    end
end
nactive=length(indices);
if(nactive>=1 && nconE>=1)
    dc_dx=dc_dx(:,indices);
    c2=c(indices);
    ConMat=[dc_dx;dceq_dx];
elseif(nactive>=1)
    dc_dx=dc_dx(:,indices);
    c2=c(indices);
    ConMat=dc_dx;
elseif(nconE>=1)
    ConMat=dceq_dx;
end
if(nactive+nconE>=1)
    lambda=-ConMat\dm_dx;
    lambda(lambda<0)=0;
    dL=dm_dx+ConMat*lambda;
else
    dL=dm_dx;
    lambda=0;
end
if(isempty(c2))
    c2=0;
end
K=max([c2;abs(ceq);0]);

s=ones(n,1);
x_hist=xk;
f_hist=fk;
d_hist=delta;
dL_hist=dL;

disp([num2str(0),': dL=',num2str(norm(dL)),', K: ',num2str(K),', fk: ',num2str(fk),', d:',num2str(delta)]);
while(~((norm(dL)<epsilon  && K<epsilon && delta<epsilon) || (norm(s)<epsilon/100 && delta<epsilon)) && (K>epsilon || fk>0) && iter<maxiter)
    % Update the finite difference length
    fd=min(epsilon/10,delta/10);
    
    % Update the penalty parameter:
    sigma=100*max(1/delta^1.1,exp(iter/10));
    
    % Update convergence tolerance:
    epsOpt=min(epsilon/100,delta/100);
    options=optimset('MaxIter',maxsubiter,'Algorithm','active-set','Display','off','TolFun',epsOpt,'TolX',epsOpt,'TolCon',epsOpt);

  
    
    % find direction sk
    
    s0=zeros(n,1);
    ConMat2=1./ConMat;
    if(K>epsilon*2 && norm(ConMat2*c2)>delta+epsOpt)
        % penalty subproblem:
        Subprob=0;

        [s,mks,exitflag,output,lambda1,dm_dx]=fmincon(@(s) PenaltyConst(xk+s,sigma,fn,Y2,Coeffs,corr,params,constr,fn_param,fid),s0,[],[],[],[],[],[],@(s) mycons(s,delta),options);
        if(norm(s)>delta+epsOpt)
            s=s*delta/norm(s);
        end
    else
        % Explicit constraint subproblem
        [s,mks,exitflag,output,lambda1,dm_dx]=fmincon(@(s) evalRBFConst(xk+s,fn,Y2,Coeffs,corr,params,fn_param,fid),s0,[],[],[],[],[],[],@(s) myconsConst(xk+s,xk,delta,constr,fn_param),options);
        Subprob=1;
        if(any(isnan(s)))
            disp('warning...warning...step is NaN')
            s=zeros(n,1);
        end
        
        if(any(isnan(s)) || norm(s)>delta+epsOpt)
            Subprob=0;
            % If that fails, revert to the penalty formulation:

            [s,mks,exitflag,output,lambda1,dm_dx]=fmincon(@(s) PenaltyConst(xk+s,sigma,fn,Y2,Coeffs,corr,params,constr,fn_param,fid),s0,[],[],[],[],[],[],@(s) mycons(s,delta),options);
            if(norm(s)>delta+epsOpt)
                s=s*delta/norm(s);
            end
        end
    end
    mks=evalRBFConst(xk+s,fn,Y2,Coeffs,corr,params,fn_param,fid);
    if(iter==0)
        s_hist=s;
    else
        s_hist=[s_hist,s];
    end
    
    % Check to see if this is a new sample point:
    i=1;
    found=-1;
    while(i<=length(points(1,:)) && found==-1)
        if(norm(points(:,i)-(xk+s))<=1e-12)
            found=i;
        end
        i=i+1;
    end
    if(found==-1)
        fks=fn(xk+s,fid(2),fn_param{:});
        flows=fn(xk+s,fid(1),fn_param{:});
    else
        flows=fn(xk+s,fid(1),fn_param{:});
        fks=flows+fvals(found);
    end
%*************************************************************
    % Compute merit function and surrogate:
    AugL=PenaltyConst(xk,sigma,fn,Y2,Coeffs,corr,params,constr,fn_param,fid);
    AugLms=PenaltyConst(xk+s,sigma,fn,Y2,Coeffs,corr,params,constr,fn_param,fid);
    AugLs=AugLms-mks+fks;  
    
%*************************************************************
    if(AugLms>=AugL-a*delta)
        rho=0;
%         disp('problem with low-fid problem')
    else
        rho=(AugL-AugLs)/(AugL-AugLms);
    end
    
    
    % update trust region radius
    if(rho>=eta2)
        delta=delta;
    elseif(rho>=eta1) 
        delta=min(gamma1*delta,deltamax);
    elseif(rho<=eta0) 
        delta=gamma0*delta;
    else
        delta=delta;
    end
    
    % update the iterate xk
    xk_old=xk;
    if(AugLs<AugL)
        xk=xk+s;
        fk=fks;
        flow=flows;
        points=[points,xk];
        fvals=[fvals;fk-flow];
        
        % Compute lagrange multipliers:
        dm_dx=zeros(n,1);
        mk=evalRBFConst(xk,fn,Y2,Coeffs,corr,params,fn_param,fid);
        [c,ceq]=constr(xk,fn_param{:});
        ncon=length(c);
        nconE=length(ceq);
        dc_dx=zeros(n,ncon);
        dceq_dx=zeros(n,nconE);
        ConMat=zeros(n,ncon+nconE);
        for i=1:n
            dx=zeros(n,1);
            dx(i)=fd;
        %     dm_dx(i)=(fn(xk+dx,fidelity)-mk)/dx(i);
            dm_dx(i)=(evalRBFConst(xk+dx,fn,Y2,Coeffs,corr,params,fn_param,fid)-mk)/(fd);
            [c2,ceq2]=constr(xk+dx,fn_param{:});
            if(ncon>=1)
                dc_dx(i,:)=(c2'-c')/fd;
            end
            if(nconE>=1)
                dceq_dx(i,:)=(ceq2'-ceq')/fd;
            end
        end

        indices=[];
        c2=0;
        for i=1:ncon
            if(c(i)>=-epsOpt && dc_dx(:,i)'*dm_dx<=0)
                indices=[indices,i];
            end
        end
        nactive=length(indices);
        if(nactive>=1 && nconE>=1)
            dc_dx=dc_dx(:,indices);
            c2=c(indices);
            ConMat=[dc_dx;dceq_dx];
        elseif(nactive>=1)
            dc_dx=dc_dx(:,indices);
            c2=c(indices);
            ConMat=dc_dx;
        elseif(nconE>=1)
            ConMat=dceq_dx;
        end
        if(nactive+nconE>=1)
            lambda=-ConMat\dm_dx;
            lambda(lambda<0)=0;
            dL=dm_dx+ConMat*lambda;
        else
            dL=dm_dx;
            lambda=0;
        end
        if(isempty(c2))
            c2=0;
        end
        K=max([c2;abs(ceq);0]);
        
        
        
    else
        xk=xk;
        fk=fk;
        points=[points,xk_old+s];
        fvals=[fvals;fks-flows];
        AugL=AugL;
%         lambda
%         dm_dx=dm_dx_old;
%         ConMat=ConMat_old;
    end

    % Update model:
    [Y2,Coeffs,linear,points,fvals,params]=RBFModelConst(fn,xk,fk-flow,points,fvals,constr,theta1,theta2,theta3,theta4,delta,deltamax,pmax,corr,params,fn_param,fid,optTheta);
    
    % Compute Lagrange multipliers:
    dm_dx=zeros(n,1);
    mk=evalRBFConst(xk,fn,Y2,Coeffs,corr,params,fn_param,fid);
    [c,ceq]=constr(xk,fn_param{:});
    ncon=length(c);
    nconE=length(ceq);
    dc_dx=zeros(n,ncon);
    dceq_dx=zeros(n,nconE);
    ConMat=zeros(n,ncon+nconE);
    for i=1:n
        dx=zeros(n,1);
        dx(i)=fd;
    %     dm_dx(i)=(fn(xk+dx,fidelity)-mk)/dx(i);
        dm_dx(i)=(evalRBFConst(xk+dx,fn,Y2,Coeffs,corr,params,fn_param,fid)-mk)/(fd);
        [c2,ceq2]=constr(xk+dx,fn_param{:});
        if(ncon>=1)
            dc_dx(i,:)=(c2'-c')/fd;
        end
        if(nconE>=1)
            dceq_dx(i,:)=(ceq2'-ceq')/fd;
        end
    end

    indices=[];
    c2=0;
    for i=1:ncon
        if(c(i)>=-epsOpt && dc_dx(:,i)'*dm_dx<=0)
            indices=[indices,i];
        end
    end
    nactive=length(indices);
    if(nactive>=1 && nconE>=1)
        dc_dx=dc_dx(:,indices);
        c2=c(indices);
        ConMat=[dc_dx;dceq_dx];
    elseif(nactive>=1)
        dc_dx=dc_dx(:,indices);
        c2=c(indices);
        ConMat=dc_dx;
    elseif(nconE>=1)
        ConMat=dceq_dx;
    end
    if(nactive+nconE>=1)
        lambda=-ConMat\dm_dx;
        lambda(lambda<0)=0;
        dL=dm_dx+ConMat*lambda;
    else
        dL=dm_dx;
        lambda=0;
    end
    if(isempty(c2))
        c2=0;
    end
    K=max([c2;abs(ceq);0]);
    iter=iter+1;
    
    if(Subprob)
        % solved explicit constraint subproblem
        disp([num2str(iter),': dL=',num2str(norm(dL)),', K: ',num2str(K),', fk: ',num2str(fk),', fks: ',num2str(fks),', d:',num2str(delta)]);
    else
        % solved penalty subproblem
        disp([num2str(iter),'*: dL=',num2str(norm(dL)),', K: ',num2str(K),', fk: ',num2str(fk),', fks: ',num2str(fks),', d:',num2str(delta)]);
    end
    x_hist=[x_hist,xk];
    f_hist=[f_hist,fk];
    d_hist=[d_hist,delta];
    dL_hist=[dL_hist,dL];
    
    % If numerical difficulties emerge with fmincon and bound minimum
    % trust-region size:
%     if(delta<epsilon/1000)
%         delta=epsilon/1000;
%     end
end


disp(' ')
disp(['Number of high-fidelity Constraint calls=',num2str(length(points)),'.'])
disp(' ');
disp(' ');
fvals=fvals';