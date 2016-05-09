% function [xk,fk,points,fvals]=unConst(fn,points,fvals,corr,params,fid,fn_param,optTheta)
%
%   This function is used to find a feasible point for the multifidelity
%   constraint optimization, essentially the algorithm in AIAA 2010-2912
%   with additional stopping criteria and the parameter "d" added to the
%   merit functions to assure the value of the constraint is bounded from
%   below.
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
function [xk,fk,points,fvals]=unConst(fn,points,fvals,corr,params,fid,fn_params,optTheta)

disp('Minimizing High-Fidelity Constraint Violation')

xk=points(:,1);
flow=fn(points(:,1)',fid(1),fn_params{:});
fk=fvals(1)+flow;
mk=fk;
n=length(xk);


fid=[2,3]; % fidelity levels contained in fn


% Algorithm Parameters:    ***********************************************
%       See AIAA-2010-2912 for details, beta and mu defined in Ref: 16
pmax=50;                    % maximum number of calibration points
eta0=0;                     % trust region movement criteria (default: 0)
eta1=.2;                    % trustregion expansion crtiteria
gamma0=1/2;                 % trust region contraction ratio
gamma1=2;                   % trust region expansion ratio
delta0=max(1,norm(xk,inf)); % Initial trust region size
deltamax=10^3*delta0;       % Maximum trust region size
epsilon=5*10^(-4);          % Convergence tolerance
kd=10^-4;                   % Fraction of Cauchy decrease constant
alpha=.9;                   % Line-search reduction ratio
mu=200;                     % termination criteria, see Ref 16
beta=3;                     % termination criteria, see Ref 16
theta1=10^-3;               % RBF Model Construction Parameter, see Ref 17
theta2=10^-4;               % RBF Model Construction Parameter, see Ref 17
theta3=10;                  % RBF Model Construction Parameter, see Ref 17
theta4=max(sqrt(n),10);     % RBF Model Construction Parameter, see Ref 17
maxiter=200;                % maximum number of iterations
maxsubiter=200;             % maximum number of iterations on TR subproblem      
fd_type=1;                  % Finite difference type, 1-first order, 2-2nd



% Algorithm:  *************************************************************

% Initial RBF Model Structures:
krig.corr=corr;
krig.params=params;
krig.points=points;
krig.Y2=[];
krig.fvals=fvals;
krig.Coeffs=[];
krig.L=[];
krig.sigma_2=0;

LowF_evals=0;

% Build RBF Error models:
delta=delta0;   
[krig,linear]=RBFModel(fn,fid,xk,[flow;fk],krig,theta1,theta2,theta3,theta4,delta,deltamax,pmax,optTheta,fn_params);

fd=min(1e-5, delta/100); % initial finite difference step,
%                       ...gets adjusted later
dm_dx=zeros(n,1);
for i=1:n
    dx=zeros(n,1);
    dx(i)=fd;
    if(fd_type==1)
        dm_dx(i)=(evalRBF(xk+dx,fn,fid,krig,fn_params)-mk)/fd;
    else
        dm_dx(i)=(evalRBF(xk+dx,fn,fid,krig,fn_params)-evalRBF(xk-dx,fn,fid,krig,fn_params))/(2*fd);
    end
end
LowF_evals=LowF_evals+n;

s=ones(n,1);
x_hist=xk;
f_hist=fk;
d_hist=delta;
dm_hist=dm_dx;
iter=0;
disp(['  0#: ck: ',num2str(fk)]);
while(~((norm(dm_dx)<epsilon || norm(s)<epsilon/100) && linear && delta<=mu*norm(dm_dx)) && fk>0 && iter<250)
    % Update Optimization parameters based on TR size:
    fd=min(1e-5, delta/100);
    epsOpt=min(epsilon/10,delta/100);
    options=optimset('MaxIter',maxsubiter,'Algorithm','active-set','Display','off','TolFun',epsOpt,'TolX',epsOpt,'TolCon',epsOpt);
    
    % Make sure trust region size is in appropriate range, see Conn Ref: 16
    % This also serves to shrink trust region to zero if model gradient is
    % approximately zero, the termination criteria in the algorithm
    if((norm(dm_dx)<epsilon || norm(s)<epsilon/100) && (~linear || delta>mu*norm(dm_dx)))

        % Run final criticality algorithm
        delta_hat=delta;
        [krig,linear]=RBFModel(fn,fid,xk,[flow;fk],krig,theta1,theta2,theta3,theta4,delta_hat,deltamax,pmax,optTheta,fn_params);
        mk=evalRBF(xk,fn,fid,krig,fn_params);
        LowF_evals=LowF_evals+1;
        dm_dx=zeros(n,1);
        for i=1:n
            dx=zeros(n,1);
            dx(i)=fd;
            if(fd_type==1)
                dm_dx(i)=(evalRBF(xk+dx,fn,fid,krig,fn_params)-mk)/fd;
            else
                dm_dx(i)=(evalRBF(xk+dx,fn,fid,krig,fn_params)-evalRBF(xk-dx,fn,fid,krig,fn_params))/(2*fd);
            end
            
        end
        LowF_evals=LowF_evals+n;
        
        % Fix size of trust region until sizing is appropriate:
        while(delta_hat>mu*norm(dm_dx) && delta_hat>=epsilon/100)  
            delta_hat=alpha*delta_hat;
            % update % update m to be linear on delta_hat
            [krig,linear]=RBFModel(fn,fid,xk,[flow;fk],krig,theta1,theta2,theta3,theta4,delta_hat,deltamax,pmax,optTheta,fn_params);
            mk=evalRBF(xk,fn,fid,krig,fn_params);
            LowF_evals=LowF_evals+1;
            dm_dx=zeros(n,1);
            for i=1:n
                dx=zeros(n,1);
                dx(i)=fd;
                if(fd_type==1)
                    dm_dx(i)=(evalRBF(xk+dx,fn,fid,krig,fn_params)-mk)/fd;
                else
                    dm_dx(i)=(evalRBF(xk+dx,fn,fid,krig,fn_params)-evalRBF(xk-dx,fn,fid,krig,fn_params))/(2*fd);
                end
            end
            LowF_evals=LowF_evals+n;
        end
        
        delta=max(delta_hat,beta*norm(dm_dx));
        delta=min(delta,deltamax);
    
    end
    
    % Find the direction sk to minimize the surrogate model:
    s0=zeros(n,1);
    [s,mks,exitflag,output,lambda,dm_dx]=fmincon(@(s) evalRBF(xk+s,fn,fid,krig,fn_params),s0,[],[],[],[],[],[],@(s) mycons(s,delta),options);
    LowF_evals=LowF_evals+output.funcCount;
    
    % Verify that fmincon hasn't died:
    if(any(isnan(s)))
        s=zeros(n,1);
    elseif(norm(s)>delta+epsOpt)
        s=s*delta/norm(s);
        mks=evalRBF(xk+s,fn,fid,krig,fn_params);
        disp('    fmincon may be having trouble');
    end
    
    % Check the fraction of Cauchy decrease condition:
    H=zeros(length(xk));
    for i=1:length(xk)
        dx=zeros(n,1);
        dx(i)=fd;
        if(fd_type==1)
            dm_dx(i)=(evalRBF(xk+dx,fn,fid,krig,fn_params)-mk)/fd;
        else
            dm_dx(i)=(evalRBF(xk+dx,fn,fid,krig,fn_params)-evalRBF(xk-dx,fn,fid,krig,fn_params))/(2*fd);
        end
        H(i,i)=(evalRBF(xk+dx,fn,fid,krig,fn_params)-2*mk+evalRBF(xk-dx,fn,fid,krig,fn_params))/fd^2;
        LowF_evals=LowF_evals+4;
        for j=1:length(xk)
            if(j~=i)
                dx2=zeros(n,1);
                dx2(j)=fd;
                % added the 2:
                dm_dx2=(evalRBF(xk+dx+dx2,fn,fid,krig,fn_params)-evalRBF(xk-dx+dx2,fn,fid,krig,fn_params))/(2*fd);
                H(i,j)=(dm_dx2-dm_dx(i))/fd;
                H(j,i)=H(i,j);
                LowF_evals=LowF_evals+2;
            end
        end
    end
    kh=norm(H)+1; % account for division by zero
    
    % Use linesearch method to satisfy cauchy decrease requirement:
    if(mk-mks<kd/2*norm(dm_dx)*min(norm(dm_dx)/kh,delta))
        LowF_evals=LowF_evals+1;
        % insufficient decrease
        disp('    insufficient decrease found, using line-search');
        s=-dm_dx/norm(dm_dx)*delta;
        while(mk-mks<kd/2*norm(dm_dx)*min(norm(dm_dx)/kh,delta) && norm(s)>=epsilon/110)
           LowF_evals=LowF_evals+1;
           s=alpha*s;
        end
        mks=evalRBF(xk+s,fn,fid,krig,fn_params);
        LowF_evals=LowF_evals+1;
    end
    if(iter==0)
        s_hist=s;
    else
        s_hist=[s_hist,s];
    end
    
    % See if xk+s has been evaluated before:
    i=1;
    found=-1;
    while(i<=length(points(1,:)) && found==-1)
        if(norm(points(:,i)-(xk+s))<=1e-12)
            found=i;
        end
        i=i+1;
    end
    % If this is a new point, evaluate it:
    if(found==-1)
        fks=fn(xk+s,fid(2),fn_params{:});
        flows=fn(xk+s,fid(1),fn_params{:});
        LowF_evals=LowF_evals+1;
    else
        flows=fn(points(:,i),fid(1),fn_params{:});
        LowF_evals=LowF_evals+1;
        fks=fmeds+krig.fvals(found);
        s=points(:,i)-xk;
    end
    
    if(mk<=mks) %Again checking to see if fmincon died
        rho=0;
    else
        rho=(fk-fks)/(mk-mks);
    end
    
    
    % update trust region radius
    if(rho>=eta1 && delta<beta*norm(dm_dx)) 
        delta=min(gamma1*delta,deltamax);
    elseif(rho>=eta1 && delta>=beta*norm(dm_dx))
        delta=delta;
    elseif(rho<eta1 && ~linear)
        delta=delta;
    else % rho<eta1 & linear
        delta=gamma0*delta;
    end
    
    % update the iterate xk
    xk_old=xk;
    if(rho>=eta1)
        xk=xk+s;
        fk=fks;
        flow=flows;
        krig.points=[krig.points,xk];
        krig.fvals=[krig.fvals;fk-flow];
    elseif(rho>eta0 && linear)
        xk=xk+s;
        fk=fks;
        flow=flows;
        krig.points=[krig.points,xk];
        krig.fvals=[krig.fvals;fk-flow];
    else
        xk=xk;
        fk=fk;
        krig.points=[krig.points,xk_old+s];
        krig.fvals=[krig.fvals;fks-flows];
    end
    
    
    
    % form new model, m_k+1
    [krig,linear]=RBFModel(fn,fid,xk,[flow;fk],krig,theta1,theta2,theta3,theta4,delta,deltamax,pmax,optTheta,fn_params);
    
    iter=iter+1;
    
    % Compute the model gradient:
    mk=evalRBF(xk,fn,fid,krig,fn_params);
    LowF_evals=LowF_evals+1;
    dm_dx=zeros(n,1);
    for i=1:n
        dx=zeros(n,1);
        dx(i)=fd;
        if(fd_type==1)
            dm_dx(i)=(evalRBF(xk+dx,fn,fid,krig,fn_params)-mk)/fd;
        else
            dm_dx(i)=(evalRBF(xk+dx,fn,fid,krig,fn_params)-evalRBF(xk-dx,fn,fid,krig,fn_params))/(2*fd);
        end
    end
    LowF_evals=LowF_evals+n;

    % Output iteration info:
    disp([num2str(iter),': fk: ',num2str(fk),' ||dm_dx||: ',num2str(norm(dm_dx)),' TR: ',num2str(delta)]);

    x_hist=[x_hist,xk];
    f_hist=[f_hist,fk];
    d_hist=[d_hist,delta];
    dm_hist=[dm_hist,dm_dx];
    
    % This has been added for numerical issues, the need really depends on
    % how if the subproblem is solved successfully:
    if(delta<epsilon/150)
        delta=epsilon/150;
    end
    
end
LowF_evals=LowF_evals+length(krig.points(1,:));

% Create these vectors for convergence plots:
s2=zeros(1,length(s_hist(1,:)));
for i=1:length(s2)
    s2(i)=norm(s_hist(:,i));
end
dm_hist2=zeros(1,length(dm_hist(1,:)));
for i=1:length(dm_hist2)
    dm_hist2(i)=norm(dm_hist(:,i));
end
dx_hist=diff(x_hist,1,2);
dx_hist2=zeros(1,length(dx_hist(1,:)));
for i=1:length(dx_hist2)
    dx_hist2(i)=norm(dx_hist(:,i));
end

disp(['   High-Fidelity Evals: ',num2str(length(krig.points(1,:))),' Low-Fidelity Evals: ',num2str(LowF_evals)])

fvals=krig.fvals';
points=krig.points;
