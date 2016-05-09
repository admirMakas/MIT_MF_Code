close all;
clc;
font=12;

% Function, Bounds, and Parameter Definition:   ***************************

fn=@analyze; % fn - a function handle to the function being optimized
constr=@myAFconst; % constr - a function handle to the function 
                   %    containing the constraints (fmincon format)

ncpts=5; % number of control points on airfoil upper or lower surface
n=2*ncpts+1; % design variables, 2*#control points+AoA
ndiv=1; % number of divisions used in initial Latin hypercube sample,
%       for none, use ndiv=1 and random starting point is selected


% Upper and Lower bounds for Latin hypercube and the optimization
LB=[-5;zeros(ncpts,1);zeros(ncpts,1)-.03];
UB=[5;zeros(ncpts,1)+.03;zeros(ncpts,1)];
scale=ones(n,1);
scale(1)=1/10;
LB=LB.*scale;
UB=UB.*scale;

% Airfoil optimization parameters
FoilGen=1; % which airfoil generator to use,
%           0 - Fourier series coefficients
%           1 - Spline points
%           2 - Points on the airfoil

M=1.5; % Mach number
gamma=1.4; % specific heat ratio
L_D_or_D=0; % Will minimize the drag coefficient if 0, L/D if 1.
fn_params={M, gamma, FoilGen,L_D_or_D,scale}; % additional parameters to fn
%           This parameters are not changed during the optimization,
%           {} for none


fid=[2,3]; % fidelity levels contained in fn



% Algorithm Parameters:    ***********************************************
%       See MDAO paper for details
pmax=20;                    % maximum number of calibration points
eta0=0.25;                  % trust region contraction criteria
eta1=.75;                   % trustregion expansion crtiteria
eta2=2;                     % trust region expansion criteria
gamma0=1/2;                 % trust region contraction ratio
gamma1=2;                   % trust region expansion ratio
delta0=max(1,norm(UB,inf)); % Initial trust region size
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

% RBF Model Parameters:   ************************************************
% Note currently only supports isometric correlation functions!
% Details are given in phi.m try: "help phi" at the command line
corr='Gauss'; % Correlation function
gammaR=2;     % gamma parameter (\xi in AIAA-2010-2912 for Gauss)
betaR=[];     % beta parameter, see phi.m
params=[betaR,gammaR];
optTheta=false;  % Boolean, true: use maximum likelihood estimator for gamma
                % false: use gamma parameter given
                % Note: this is done by checking 10-20 values and using the
                % one with the highest likelihood, see AddPoints.m

% Initial latin hypercube sample or random starting point for algorithm****
LHS=lhsdesign(n,ndiv);
points=zeros(n,ndiv);
f1=zeros(ndiv,1);
f2=zeros(ndiv,1);
for i=1:ndiv
    points(:,i)=LB+(UB-LB).*LHS(:,i);
    f1(i)=fn(points(:,i)',fid(1),fn_param{:});
    f2(i)=fn(points(:,i)',fid(2),fn_param{:});
end
fvals=f2-f1;

i=1;
Xo=points(:,i);

LB=[];
UB=[];
if(i<=ndiv)
    xk=points(:,i);
    fk=f2(i);
    flow=f1(i);
    mk=f2(i);
else
    xk=points(:,1);
    fk=f2(1);
    flow=f1(1);
    mk=f2(1);
end
iter=0;
delta=delta0;
epsOpt=epsilon/10;

% Build initial surrogate model for the objective function:
[Y2,Coeffs,linear,points,fvals,params]=RBFModel1(fn,xk,fk-flow,points,fvals,constr,theta1,theta2,theta3,theta4,delta,deltamax,pmax,corr,params,fn_param,fid,optTheta);


% Compute the initial Lagrange multipliers:
dm_dx=zeros(n,1);
mk=evalRBF(xk,fn,Y2,Coeffs,corr,params,fn_param,fid);
[c,ceq]=constr(xk,fn_param{:});
ncon=length(c);
nconE=length(ceq);
dc_dx=zeros(n,ncon);
dceq_dx=zeros(n,nconE);
ConMat=zeros(n,ncon+nconE);
for i=1:n
    dx=zeros(n,1);
    dx(i)=fd;
    dm_dx(i)=(evalRBF(xk+dx,fn,Y2,Coeffs,corr,params,fn_param,fid)-mk)/(fd);
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
while(~((norm(dL)<epsilon  && K<epsilon && delta<epsilon) || (norm(s)<epsilon/100 && delta<epsilon)) && iter<maxiter)
    % Update the finite difference length
    fd=min(epsilon/10,delta/10);
    
    % Update the penalty parameter:
    sigma=max(1/delta^1.1,exp(iter/10))+100;
    
    % Update convergence tolerance:
    epsOpt=min(epsilon/100,delta/100);
    options=optimset('MaxIter',maxsubiter,'Algorithm','active-set','Display','off','TolFun',epsOpt,'TolX',epsOpt,'TolCon',epsOpt);

  
    
    % find direction sk
    
    s0=zeros(n,1);
    ConMat2=1./ConMat;
    if(K>epsilon*2 && norm(ConMat2*c2)>delta+epsOpt)
        % penalty subproblem:
        Subprob=0;
%         [s,mks,exitflag,output,lambda1,dm_dx]=fmincon(@(s) Penalty(xk+s,sigma,fn,Y2,Coeffs,corr,params,constr,fn_param,fid),s0,[],[],[],[],[],[],@(s) mycons(xk+s,xk,delta,constr,fn_param),options);
        [s,mks,exitflag,output,lambda1,dm_dx]=fmincon(@(s) Penalty(xk+s,sigma,fn,Y2,Coeffs,corr,params,constr,fn_param,fid),s0,[],[],[],[],[],[],@(s) TR(xk+s,xk,delta),options);
        if(norm(s)>delta+epsOpt)
            s=s*delta/norm(s);
        end
    else
        % Explicit constraint subproblem
        [s,mks,exitflag,output,lambda1,dm_dx]=fmincon(@(s) evalRBF(xk+s,fn,Y2,Coeffs,corr,params,fn_param,fid),s0,[],[],[],[],[],[],@(s) mycons(xk+s,xk,delta,constr,fn_param),options);
        Subprob=1;
        if(any(isnan(s)))
            disp('warning...warning...step is NaN')
            s=zeros(n,1);
        end
        
        if(any(isnan(s)) || norm(s)>delta+epsOpt)
            Subprob=0;
            % If that fails, revert to the penalty formulation:
%             [s,mks,exitflag,output,lambda1,dm_dx]=fmincon(@(s) Penalty(xk+s,sigma,fn,Y2,Coeffs,corr,params,constr,fn_param,fid),s0,[],[],[],[],[],[],@(s) mycons(xk+s,xk,delta,constr,fn_param),options);
            [s,mks,exitflag,output,lambda1,dm_dx]=fmincon(@(s) Penalty(xk+s,sigma,fn,Y2,Coeffs,corr,params,constr,fn_param,fid),s0,[],[],[],[],[],[],@(s) TR(xk+s,xk,delta),options);
            if(norm(s)>delta+epsOpt)
                s=s*delta/norm(s);
            end
        end
    end
    mks=evalRBF(xk+s,fn,Y2,Coeffs,corr,params,fn_param,fid);
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
    AugL=Penalty(xk,sigma,fn,Y2,Coeffs,corr,params,constr,fn_param,fid);
    AugLms=Penalty(xk+s,sigma,fn,Y2,Coeffs,corr,params,constr,fn_param,fid);
    AugLs=AugLms-mks+fks;  
    
%*************************************************************
    if(AugLms>=AugL-a*delta)
        rho=0;
%         disp('problem with low-fid problem')
    else
        rho=(AugL-AugLs)/(AugL-AugLms);
    end
    
    
    % update trust region radius
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
        mk=evalRBF(xk,fn,Y2,Coeffs,corr,params,fn_param,fid);
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
            dm_dx(i)=(evalRBF(xk+dx,fn,Y2,Coeffs,corr,params,fn_param,fid)-mk)/(fd);
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
    [Y2,Coeffs,linear,points,fvals,params]=RBFModel1(fn,xk,fk-flow,points,fvals,constr,theta1,theta2,theta3,theta4,delta,deltamax,pmax,corr,params,fn_param,fid,optTheta);
    
    % Compute Lagrange multipliers:
    dm_dx=zeros(n,1);
    mk=evalRBF(xk,fn,Y2,Coeffs,corr,params,fn_param,fid);
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
        dm_dx(i)=(evalRBF(xk+dx,fn,Y2,Coeffs,corr,params,fn_param,fid)-mk)/(fd);
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
s2=zeros(1,length(s_hist(1,:)));
for i=1:length(s2)
    s2(i)=norm(s_hist(:,i));
end
dL_hist2=zeros(1,length(dL_hist(1,:)));
for i=1:length(dL_hist2)
    dL_hist2(i)=norm(dL_hist(:,i));
end
dx_hist=diff(x_hist,1,2);
dx_hist2=zeros(1,length(dx_hist(1,:)));
for i=1:length(dx_hist2)
    dx_hist2(i)=norm(dx_hist(:,i));
end


% Plot the convergence rates:
figure
h=gca;
set(h,'FontSize',font)
semilogy(abs(f_hist))
hold on;
semilogy(abs(f_hist-f_hist(end)),'r')
semilogy(abs(d_hist),'g')
% semilogy(s2,':k')
semilogy(dL_hist2,'k')
% semilogy(dx_hist2,'c')
hold off;
a=get(h,'Children');
set(a,'LineWidth',2)
xlabel('Iteration')
ylabel('Convergence Quantity')
title(['Convergence History, ',num2str(length(points(1,:))),' Function Calls'])
legend('f_{high}(x_i)','|f_{high}(x_i)-f_{high}(x^*)|','Trust-Region Radius','||dL/dx||',3)

disp(['Number of high-fidelity calls=',num2str(length(points)),'.'])



% Plot the airfoil results:
if(n>=3)
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
    title(['Optimized Airfoil, D=',num2str(f_hist(end)),' \alpha=',num2str(xk(1)*scale(1)),'^\circ'])
    xlabel('x/c')
    ylabel('t/c')
    a=get(h,'Children');
    set(a,'LineWidth',2,'MarkerSize',5)
end