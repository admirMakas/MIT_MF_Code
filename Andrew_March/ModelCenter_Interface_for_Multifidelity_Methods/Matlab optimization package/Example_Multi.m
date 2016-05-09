close all;
clc
font=12;


% Function, Bounds, and Parameter Definition:   ***************************

% fn - a function handle to the function being optimized
fn=@ModelCenterWrap2;


ndiv=1; % number of divisions used in initial Latin hypercube sample,
%       for none, use ndiv=1 and random starting point is selected
n=2;

% Upper and Lower bounds for Latin hypercube and the optimization
LB=[-5;-5];
UB=[5;5];
scale=ones(n,1);
LB=LB.*scale;
UB=UB.*scale;


fn_params={}; % additional parameters to fn
%           This parameters are not changed during the optimization,
%           {} for none


fid=[0,1,2]; % fidelity levels contained in fn


% Algorithm Parameters:    ***********************************************
%       See AIAA-2010-2912 for details, beta and mu defined in Ref: 16
pmax=50;                    % maximum number of calibration points
eta0=0;                     % trust region movement criteria (default: 0)
eta1=.2;                    % trustregion expansion crtiteria
gamma0=1/2;                 % trust region contraction ratio
gamma1=2;                   % trust region expansion ratio
delta0=max(1,norm(UB,inf)); % Initial trust region size
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

% RBF Model Parameters:   ************************************************
% Note currently only supports isometric correlation functions!
% Details are given in phi.m try: "help phi" at the command line
corr='Gauss'; % Correlation function
gammaR=2;     % gamma parameter (\xi in AIAA-2010-2912 for Gauss)
betaR=[];     % beta parameter, see phi.m
params=[betaR,gammaR];
optTheta=true;  % Boolean, true: use maximum likelihood estimator for gamma
                % false: use gamma parameter given
                % Note: this is done by checking 10-20 values and using the
                % one with the highest likelihood, see AddPoints.m


% Initial latin hypercube sample or random starting point for algorithm****
LowF_evals=0;
LHS=lhsdesign(n,ndiv);
points=zeros(n,ndiv);
f1=zeros(ndiv,1);
f2=zeros(ndiv,1);
f3=zeros(ndiv,1);
for i=1:ndiv
    points(:,i)=LB+(UB-LB).*LHS(:,i);
    f1(i)=fn(points(:,i),fid(1),fn_params{:});
    f2(i)=fn(points(:,i),fid(2),fn_params{:});
    f3(i)=fn(points(:,i),fid(3),fn_params{:});
end
fvalsHM=f3-f2;
fvalsHL=f3-f1;
LowF_evals=LowF_evals+ndiv; % counter for low-fidelity evaluations


% Algorithm:  *************************************************************

% Initial RBF Model Structures:
krigM.corr=corr;
krigM.params=params;
krigM.points=points;
krigM.Y2=[];
krigM.fvals=fvalsHM;
krigM.Coeffs=[];
krigM.L=[];
krigM.sigma_2=0;

krigL.corr=corr;
krigL.params=params;
krigL.points=points;
krigL.Y2=[];
krigL.fvals=fvalsHL;
krigL.Coeffs=[];
krigL.L=[];
krigL.sigma_2=0;

xk=points(:,1);
fk=f3(1);
flow=f1(1);
fmed=f2(1);
mk=f3(1);

% Build RBF Error models:
delta=delta0;   
[krigM,krigL,linear]=RBFModelm(fn,fid,xk,[flow;fmed;fk],krigM,krigL,theta1,theta2,theta3,theta4,delta,deltamax,pmax,optTheta,fn_params);

fd=min(1e-5, delta/100); % initial finite difference step,
%                       ...gets adjusted later
dm_dx=zeros(n,1);
for i=1:n
    dx=zeros(n,1);
    dx(i)=fd;
    if(fd_type==1)
        dm_dx(i)=(evalRBFm(xk+dx,fn,fid,krigM,krigL,fn_params)-mk)/fd;
    else
        dm_dx(i)=(evalRBFm(xk+dx,fn,fid,krigM,krigL,fn_params)-evalRBFm(xk-dx,fn,fid,krigM,krigL,fn_params))/(2*fd);
    end
end
LowF_evals=LowF_evals+n;

s=ones(n,1);
x_hist=xk;
f_hist=fk;
d_hist=delta;
dm_hist=dm_dx;
iter=0;
while(~((norm(dm_dx)<epsilon || (norm(s)<epsilon/100 && delta<=epsilon)) && linear && delta<=mu*norm(dm_dx)) && iter<maxiter)
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
        [krigM,krigL,linear]=RBFModelm(fn,fid,xk,[flow;fmed;fk],krigM,krigL,theta1,theta2,theta3,theta4,delta_hat,deltamax,pmax,optTheta,fn_params);
        mk=evalRBFm(xk,fn,fid,krigM,krigL,fn_params);
        LowF_evals=LowF_evals+1;
        dm_dx=zeros(n,1);
        for i=1:n
            dx=zeros(n,1);
            dx(i)=fd;
            if(fd_type==1)
                dm_dx(i)=(evalRBFm(xk+dx,fn,fid,krigM,krigL,fn_params)-mk)/fd;
            else
                dm_dx(i)=(evalRBFm(xk+dx,fn,fid,krigM,krigL,fn_params)-evalRBFm(xk-dx,fn,fid,krigM,krigL,fn_params))/(2*fd);
            end
            
        end
        LowF_evals=LowF_evals+n;
        
        % Fix size of trust region until sizing is appropriate:
        while(delta_hat>mu*norm(dm_dx) && delta_hat>=epsilon/100)  
            delta_hat=alpha*delta_hat;
            % update % update m to be linear on delta_hat
            [krigM,krigL,linear]=RBFModelm(fn,fid,xk,[flow;fmed;fk],krigM,krigL,theta1,theta2,theta3,theta4,delta_hat,deltamax,pmax,optTheta,fn_params);
            mk=evalRBFm(xk,fn,fid,krigM,krigL,fn_params);
            LowF_evals=LowF_evals+1;
            dm_dx=zeros(n,1);
            for i=1:n
                dx=zeros(n,1);
                dx(i)=fd;
                if(fd_type==1)
                    dm_dx(i)=(evalRBFm(xk+dx,fn,fid,krigM,krigL,fn_params)-mk)/fd;
                else
                    dm_dx(i)=(evalRBFm(xk+dx,fn,fid,krigM,krigL,fn_params)-evalRBFm(xk-dx,fn,fid,krigM,krigL,fn_params))/(2*fd);
                end
            end
            LowF_evals=LowF_evals+n;
        end
        
        delta=max(delta_hat,beta*norm(dm_dx));
        delta=min(delta,deltamax);
    
    end
    
    % Find the direction sk to minimize the surrogate model:
    s0=zeros(n,1);
    [s,mks,exitflag,output,lambda,dm_dx]=fmincon(@(s) evalRBFm(xk+s,fn,fid,krigM,krigL,fn_params),s0,[],[],[],[],[],[],@(s) mycons(s,delta),options);
    LowF_evals=LowF_evals+output.funcCount;
    
    % Verify that fmincon hasn't died:
    if(any(isnan(s)))
        s=zeros(n,1);
    elseif(norm(s)>delta+epsOpt)
        s=s*delta/norm(s);
        mks=evalRBFm(xk+s,fn,fid,krigM,krigL,fn_params);
        disp('    fmincon may be having trouble');
    end
    
    % Check the fraction of Cauchy decrease condition:
    H=zeros(length(xk));
    for i=1:length(xk)
        dx=zeros(n,1);
        dx(i)=fd;
        if(fd_type==1)
            dm_dx(i)=(evalRBFm(xk+dx,fn,fid,krigM,krigL,fn_params)-mk)/fd;
        else
            dm_dx(i)=(evalRBFm(xk+dx,fn,fid,krigM,krigL,fn_params)-evalRBFm(xk-dx,fn,fid,krigM,krigL,fn_params))/(2*fd);
        end
        H(i,i)=(evalRBFm(xk+dx,fn,fid,krigM,krigL,fn_params)-2*mk+evalRBFm(xk-dx,fn,fid,krigM,krigL,fn_params))/fd^2;
        LowF_evals=LowF_evals+4;
        for j=1:length(xk)
            if(j~=i)
                dx2=zeros(n,1);
                dx2(j)=fd;
                % added the 2:
                dm_dx2=(evalRBFm(xk+dx+dx2,fn,fid,krigM,krigL,fn_params)-evalRBFm(xk-dx+dx2,fn,fid,krigM,krigL,fn_params))/(2*fd);
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
        mks=evalRBFm(xk+s,fn,fid,krigM,krigL,fn_params);
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
        fks=fn(xk+s,fid(3),fn_params{:});
        flows=fn(xk+s,fid(1),fn_params{:});
        fmeds=fn(xk+s,fid(2),fn_params{:});
        LowF_evals=LowF_evals+1;
    else
        flows=fn(points(:,i),fid(1),fn_params{:});
        fmeds=fn(points(:,i),fid(2),fn_params{:});
        LowF_evals=LowF_evals+1;
        fks=fmeds+krigM.fvals(found);
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
        fmed=fmeds;
        krigM.points=[krigM.points,xk];
        krigM.fvals=[krigM.fvals;fk-fmed];
        krigL.points=[krigL.points,xk];
        krigL.fvals=[krigL.fvals;fk-flow];
    elseif(rho>eta0 && linear)
        xk=xk+s;
        fk=fks;
        flow=flows;
        fmed=fmeds;
        krigM.points=[krigM.points,xk];
        krigM.fvals=[krigM.fvals;fk-fmed];
        krigL.points=[krigL.points,xk];
        krigL.fvals=[krigL.fvals;fk-flow];
    else
        xk=xk;
        fk=fk;
        krigM.points=[krigM.points,xk_old+s];
        krigM.fvals=[krigM.fvals;fks-fmeds];
        krigL.points=[krigL.points,xk_old+s];
        krigL.fvals=[krigL.fvals;fks-flows];
    end
    
    
    
    % form new model, m_k+1
    [krigM,krigL,linear]=RBFModelm(fn,fid,xk,[flow;fmed;fk],krigM,krigL,theta1,theta2,theta3,theta4,delta,deltamax,pmax,optTheta,fn_params);
    
    iter=iter+1;
    
    % Compute the model gradient:
    mk=evalRBFm(xk,fn,fid,krigM,krigL,fn_params);
    LowF_evals=LowF_evals+1;
    dm_dx=zeros(n,1);
    for i=1:n
        dx=zeros(n,1);
        dx(i)=fd;
        if(fd_type==1)
            dm_dx(i)=(evalRBFm(xk+dx,fn,fid,krigM,krigL,fn_params)-mk)/fd;
        else
            dm_dx(i)=(evalRBFm(xk+dx,fn,fid,krigM,krigL,fn_params)-evalRBFm(xk-dx,fn,fid,krigM,krigL,fn_params))/(2*fd);
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
LowF_evals=LowF_evals+length(krigM.points(1,:));

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

disp(['High-Fidelity Evals: ',num2str(length(krigM.points(1,:))),' Low-Fidelity Evals: ',num2str(LowF_evals)])


% Plot Convergence Info:  **********************************************
h=gca;
set(h,'FontSize',font)
semilogy(abs(f_hist));
hold on;
semilogy(abs(f_hist-f_hist(end)),'r')
% semilogy(s2,':k')
semilogy(dm_hist2,'k')
semilogy(abs(d_hist),'g')
hold off;
a=get(h,'Children');
set(a,'LineWidth',2)
xlabel('Iteration')
ylabel('Convergence Quantity')
title('Convergence History')
legend('|f_{high}(x_i)|','|f_{high}(x_i)-f_{high}(x*)|','||dm/dx||','Trust-region size')
a=get(h,'Children');
set(a,'LineWidth',2,'MarkerSize',5)


% Plot Initial and Final Airfoils:  ***************************************
switch(FoilGen)
    case 1
        [xinit,top0,bot0]=GenFoilSpline(points(:,1)',false);
        [xf,top,bot]=GenFoilSpline(xk',false);
    case 2
        [xinit,top0,bot0]=GenFoilPts(points(:,1)',false);
        [xf,top,bot]=GenFoilPts(xk',false);
end
 
ntop=floor(length(xk)/2);
nbot=ntop;
xctop=linspace(0,1,ntop+2)';
xcbot=linspace(0,1,nbot+2)';
% Use either...
Tpts=[0,xk(2:1+ntop)',0];
Bpts=[0,xk(2+ntop:ntop+nbot+1)',0];
Tpts0=[0,points(2:1+ntop,1)',0];
Bpts0=[0,points(2+ntop:ntop+nbot+1,1)',0];

figure
subplot(1,2,1)
h=gca;
set(h,'FontSize',font);
plot(xinit,top0,'r');
hold on;
plot(xinit,bot0,'r');
plot(xctop,Tpts0,'ob');
plot(xcbot,Bpts0,'ob');
hold off;
title('Initial Airfoil')
xlabel('x/c')
ylabel('t/c')
a=get(h,'Children');
set(a,'LineWidth',2,'MarkerSize',5)

subplot(1,2,2)
h=gca;
set(h,'FontSize',font);
plot(xf,top,'r');
hold on;
plot(xf,bot,'r');
plot(xctop,Tpts,'ob');
plot(xcbot,Bpts,'ob');
hold off;
title(['Optimized Airfoil, D=',num2str(f_hist(end)),' \alpha=',num2str(xk(1)*10),'^\circ'])
xlabel('x/c')
ylabel('t/c')
a=get(h,'Children');
set(a,'LineWidth',2,'MarkerSize',5)