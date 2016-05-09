close all;
clc;
font=12;

% Function, Bounds, and Parameter Definition:   ***************************
fn=@LDest;      % fn - a function handle to the function being optimized
constr=@Dmax;   % constr - a function handle to the function 
                   %    containing the constraint that will be approximated
constrNoAppr=@myAFconst; % Function handle to the function constainting the
                         %   containing the constraints that will not be
                         %   approximated, functions g(x) and h(x).
                         %   [] - for none

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


fid=[2,3]; % fidelity levels contained in constr.

drag=0.01;  % the value of drag in the multifidelity constraint
fn_param={M,gamma,FoilGen,ncpts,scale}; % additional parameters to fn
%           This parameters are not changed during the optimization,
%           {} for none
constr_param={drag,M,gamma,FoilGen,ncpts,scale}; % The same but for the 
%            constraint.

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
epsilon2=epsilon;           % Convergence tolerance for trust region size
a=10^-4;                    % Sufficient decrease condition
alpha=.9;                   % Line-search reduction ratio
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


                
% Beginning of Method: **********************************************
otherConstrs=isa(constrNoAppr, 'function_handle');

count=0;
countf=0;
LHS=lhsdesign(n,ndiv);
points=zeros(n,ndiv);
f1=zeros(ndiv,1);
f2=zeros(ndiv,1);
for i=1:ndiv
    points(:,i)=LB+(UB-LB).*LHS(:,i);
    if(i==1)
        ctemp=constr(points(:,i),fid(1),constr_param{:});
        nconstrs=length(ctemp);
        c1=zeros(nconstrs,ndiv);
        c2=zeros(nconstrs,ndiv);
        c1(:,i)=ctemp;
    else
        c1(:,i)=constr(points(:,i),fid(1),constr_param{:});
    end
    c2(:,i)=constr(points(:,i),fid(2),constr_param{:});  
    f1(i)=fn(points(:,i),fid(1),fn_param{:});
    f2(i)=fn(points(:,i),fid(2),fn_param{:});
end
pointsf=points;

Xo=points(:,1);
ceval=c2(:,1);
cvals=c2-c1;
fvals=f2-f1;
xk=Xo;
fk=f2(1); 
flow=f1(1);

if(otherConstrs)
    [g,h]=constrNoAppr(xk,constr_param{:});
    if(~isempty(g) && any(g>0))
        g=g(g>0);
    else
        g=0;
    end
    if(isempty(h))
        h=0;
    end
    % add inequality and equality constraints that aren't approx
else
    h=0;
    g=0;
end


iter=0;
epsOpt=epsilon/100;
epsCon=epsilon/100;
delta=delta0;
xk=Xo;
x_hist=xk;
f_hist=fk;
d_hist=delta;

K=max([g;abs(h);max(ceval,0)]);
K_hist=K;

% Find a high-fidelity constraint feasibile point:
if(ceval>0 && ~otherConstrs)
    [xk,c2,points,cvals]=unConst(constr,points,cvals,corr,params,fid,constr_param,optTheta);
    fk=fn(xk,fid(2),fn_param{:});
    flow=fn(xk,fid(1),fn_param{:});
    fvals=[fvals;fk-flow];
    pointsf=[pointsf,xk];
elseif(K>epsilon)
    [xk,c2,points,cvals]=ConstSolver(constr,points,cvals,corr,params,fid,constr_param,constrNoAppr,optTheta);
    [g,h]=constrNoAppr(xk,constr_param{:});
    if(~isempty(g) && any(g>0))
        g=g(g>0);
    else
        g=0;
    end
    if(isempty(h))
        h=0;
    end
    fk=fn(xk,fid(2),fn_param{:});
    flow=fn(xk,fid(1),fn_param{:});
    fvals=[fvals;fk-flow];
    pointsf=[pointsf,xk];
end
ceval=c2;
disp('Beginning Constrained Optimization')

% Create Initial Kriging models:
clear krigs;
krigs(nconstrs,1)= struct('Y2',[],'Coeffs',[],'params',[betaR,gammaR],'linear',true,'corr',corr);
for i=1:nconstrs
    [krig,points,cvals]=RBFModel1(constr,xk,c2-c1,points,cvals,i,theta1,theta2,theta3,theta4,delta,deltamax,pmax,corr,params,constr_param,fid,optTheta);
    krigs(i)=krig;
end

[Y2,Coeffs,linear,pointsf,fvals,params]=RBFModelConst(fn,xk,fk-flow,pointsf,fvals,constr,theta1,theta2,theta3,theta4,delta,deltamax,pmax,corr,params,fn_param,fid,optTheta);


s=ones(n,1);
LB=[];
UB=[];
lambda1=10;
dm_dx=zeros(n,1);
[mineq,meq]=SurrCons(xk,s,-1,constr,otherConstrs,constrNoAppr,krigs,constr_param,fid);

% Compute surrogate Lagrange multipliers:
if(~isempty(mineq))
    Dineq_dx=zeros(length(mineq),n);
end

if(~isempty(meq))
    Deq_dx=zeros(length(meq),n);
end
mk=fk;
for i=1:n
    dx=zeros(n,1);
    dx(i)=fd;
    % change to m!!!
    dm_dx(i)=(evalRBFFunc(xk+dx,fn,Y2,Coeffs,corr,params,fn_param,fid)-mk)/(fd);
    [mineq_p,meq_p]=SurrCons(xk+dx,s,-1,constr,otherConstrs,constrNoAppr,krigs,constr_param,fid);
    if(~isempty(mineq))
        Dineq_dx(:,i)=(mineq_p-mineq)./fd;
    end
    if(~isempty(meq))
        Deq_dx(:,i)=(meq_p-meq)./fd;
    end
end
if(~isempty(mineq))
    indices=[];
    for i=1:length(mineq)
        if(mineq(i)>=-epsOpt && Dineq_dx(i,:)*dm_dx<=0)
            indices=[indices,i];
        end
    end
    if(isempty(meq))
        A=Dineq_dx(indices,:);
    else
        A=[Dineq_dx(indices,:);Deq_dx];
    end
    A2=1./A';
    step=A2*[mineq(indices);meq]/n;
else % check for existence
    A=Deq_dx;
    A2=1./A';
    step=A2*meq/n;
end
lambda=-(A*A')\(A*dm_dx);
L=norm(dm_dx+A'*lambda);
A2=1./A';
dL_hist=L;

K=max([g;abs(h);max(ceval,0)]);

if(K<=epsilon)
    disp([num2str(iter),': dL=',num2str(L),', f: ',num2str(fk),', c: ',num2str(K),', delta=',num2str(delta)]);
else
    disp([num2str(iter),'*: dL=',num2str(L),', f: ',num2str(fk),', c: ',num2str(K),', delta=',num2str(delta)]);
end


% Begin iterations:
while(~((L<=epsilon && delta<=epsilon && K<=epsilon/10 && lambdaTR==0) || (delta<=epsilon2 && norm(s)<epsilon2/10 && K<=epsilon))  && iter<maxiter)
    % Update Finite difference
    fd=min(epsilon/10,delta/10);
    
    % Update penalty parameter
    sigma=1*max(1/delta^1.1,exp(iter/10));

    % Update convergence tolerances:
    epsOpt=min(epsilon/100,delta/100);
    epsCon=min(epsilon/100,delta/100);
    options1=optimset('MaxIter',maxsubiter,'MaxFunEval',3000,'Algorithm','active-set','Display','iter','TolFun',epsOpt,'TolX',epsOpt,'TolCon',epsOpt);
    options2=optimset('MaxIter',maxsubiter,'MaxFunEval',3000,'Algorithm','active-set','Display','iter','TolFun',epsOpt,'TolX',epsOpt,'TolCon',epsOpt);
    
    
    % find direction sk
    s0=zeros(n,1);
    [mineq,meq]=SurrCons(xk,s0,delta,constr,otherConstrs,constrNoAppr,krigs,constr_param,fid);
    Khat=max([mineq;abs(meq)]);

    [s,fks,flag,output,lambdaFmin,grad]=fmincon(@(s) evalRBFFunc(xk+s,fn,Y2,Coeffs,corr,params,fn_param,fid),s0,[],[],[],[],LB,UB,@(s) SurrCons(xk,s,delta,constr,otherConstrs,constrNoAppr,krigs,constr_param,fid),options2);

    lambda1=lambdaFmin.ineqnonlin;
    lambdaTR=lambdaFmin.ineqnonlin(end);
    
    % Verify fmincon solution is good:
    if(any(isnan(s)))
        s=zeros(n,1);
    end
    if(norm(s)>=delta+epsCon )
        disp(['problem: s>=TR size, ||s||:',num2str(norm(s)),', delta: ',num2str(delta),',Constr Vio ',output.constrviolation]);
        s=s*delta/norm(s);
    end    
    
    % store results in history:
    if(iter==0)
        s_hist=s;
        lambda_hist=lambda1(1);
    else
        s_hist=[s_hist,s];
        lambda_hist=[lambda_hist,lambda1(1)];
    end
    
%*************************************************************

    % Check if xk+sk is a new point, if not, then evaluate it:
    % For constraints:
    i=1;
    found=-1;
    while(i<=length(points(1,:)) && found==-1)
        if(norm(points(:,i)-(xk+s))<=1e-10)
            found=i;
        end
        i=i+1;
    end
    if(found==-1)
        c1s=constr(xk+s,fid(1),constr_param{:});
        c2s=constr(xk+s,fid(2),constr_param{:});
        points=[points,xk+s];
        cvals=[cvals,c2s-c1s];
    else
        c1s=constr(points(:,found),fid(1),constr_param{:});
        c2s=c1s+cvals(found);
    end
    
    % For objective function:
    i=1;
    found=-1;
    while(i<=length(pointsf(1,:)) && found==-1)
        if(norm(pointsf(:,i)-(xk+s))<=1e-10)
            found=i;
        end
        i=i+1;
    end
    if(found==-1)
        flows=fn(xk+s,fid(1),fn_param{:});
        fks=fn(xk+s,fid(2),fn_param{:});
        pointsf=[pointsf,xk+s];
        fvals=[fvals;fks-flows];
    else
        flows=fn(pointsf(:,found),fid(1),fn_param{:});
        fks=flows+fvals(found);
    end
    
    if(otherConstrs)
        [gs,hs]=constrNoAppr(xk+s,constr_param{:});
        if(~isempty(gs) && any(gs>0))
            gs=gs(gs>0);
        else
            gs=0;
        end
        if(isempty(hs))
            hs=0;
        end
        % add inequality and equality constraints that aren't approx
    else
        hs=0;
        gs=0;
    end
     
    % Verify still feasible and sufficient decrease condition:
    if(c2s<=0 && fk-fks>=a*delta) %% will need to account for approx function here..
        rho=1;
    else
        rho=0;
    end
    
    disp(['c2: ',num2str(c2),' c2s: ',num2str(c2s),' fk: ',num2str(fk),' fks: ',num2str(fks),' gs: ',num2str(max(gs))]);
    
    xk_old=xk;
    % Update if this is an improving point:
    if(fks<fk && c2s<=0  && ~any(hs>epsilon) && ~any(gs>epsilon))
        xk=xk+s;
        fk=fks;
        flow=flows;
        c2=c2s;
        c1=c1s;
        g=gs;
        h=hs;
        ceval=c2s;
        K=max([g;abs(h);max(ceval,0)]);
        
    end
    

% Update the trust region size
    if(rho>=eta1 && rho<=eta2)
        delta=min(gamma1*delta,deltamax);
    elseif(rho<=eta0)
        delta=gamma0*delta;
    else
        delta=delta;
    end

    
    % Update the Kriging models for the constraints and objectives:
    for i=1:nconstrs
        [krig,points,cvals]=RBFModel1(constr,xk,c2-c1,points,cvals,i,theta1,theta2,theta3,theta4,delta,deltamax,pmax,corr,params,constr_param,fid,optTheta);
        krigs(i)=krig;
    end
    [Y2,Coeffs,linear,pointsf,fvals,params]=RBFModelConst(fn,xk,fk-flow,pointsf,fvals,constr,theta1,theta2,theta3,theta4,delta,deltamax,pmax,corr,params,fn_param,fid,optTheta);
    
    iter=iter+1;
    
    % Compute surrogate Lagrange multipliers:
    mk=fk;
    dm_dx=zeros(n,1);
    [mineq,meq]=SurrCons(xk,s,-1,constr,otherConstrs,constrNoAppr,krigs,constr_param,fid);
    if(~isempty(mineq))
        Dineq_dx=zeros(length(mineq),n);
    end
    if(~isempty(meq))
        Deq_dx=zeros(length(meq),n);
    end
    for i=1:n
        dx=zeros(n,1);
        dx(i)=fd;
        dm_dx(i)=(evalRBFFunc(xk+dx,fn,Y2,Coeffs,corr,params,fn_param,fid)-mk)/(fd);
        [mineq_p,meq_p]=SurrCons(xk+dx,s,-1,constr,otherConstrs,constrNoAppr,krigs,constr_param,fid);
        if(~isempty(mineq))
            Dineq_dx(:,i)=(mineq_p-mineq)./fd;
        end
        if(~isempty(meq))
            Deq_dx(:,i)=(meq_p-meq)./fd;
        end
    end
    if(~isempty(mineq))
        indices=[];
        for i=1:length(mineq)
            if(mineq(i)>=-epsOpt && Dineq_dx(i,:)*dm_dx<=0)
                indices=[indices,i];
            end
        end
        if(isempty(meq))
            A=Dineq_dx(indices,:);
        else
            A=[Dineq_dx(indices,:);Deq_dx];
        end
        A2=1./A';
        step=A2*[mineq(indices);meq]/n;
    else % check for empty...
        A=Deq_dx;
        A2=1./A';
        step=A2*meq/n;
    end
    lambda=-(A*A')\(A*dm_dx);
    
    L=norm(dm_dx+A'*lambda);
    
    
    
    % Display iterations:
    disp([num2str(iter),': dL=',num2str(L),', f: ',num2str(fk),', c: ',num2str(K),', delta=',num2str(delta)]);
    
    x_hist=[x_hist,xk];
    f_hist=[f_hist,fk];
    d_hist=[d_hist,delta];
    dL_hist=[dL_hist,L];
    K_hist=[K_hist,K];

end

 % Output number of function calls:
disp(['Number of high-fidelity Objective calls=',num2str(length(fvals)),'.'])
disp(['Number of high-fidelity Constraint calls=',num2str(length(cvals(1,:))),'.'])
lambda1=lambda1(1:length(lambda1)-1);


% Make convergence history plots:***************************************
s2=zeros(1,length(s_hist(1,:)));
for i=1:length(s2)
    s2(i)=norm(s_hist(:,i));
end
x_hist2=zeros(1,length(x_hist(1,:)));
for i=1:length(x_hist)
    x_hist2(i)=norm(x_hist(:,i)-x_hist(:,end));
end
dx_hist=diff(x_hist,1,2);
dx_hist2=zeros(1,length(dx_hist(1,:)));
for i=1:length(dx_hist2)
    dx_hist2(i)=norm(dx_hist(:,i));
end

figure
h=gca;
set(h,'FontSize',font)
semilogy(abs(f_hist))
hold on;
semilogy(K_hist,'k')
semilogy(dL_hist,'r');
semilogy(x_hist2,'c')
% semilogy(dx_hist2,'r')
semilogy(abs(d_hist),'g')
% semilogy(s2,':k')
% semilogy(dx_hist2,'c')
%semilogy(abs(lambda_hist-lambda_hist(end)),'m')
hold off;
a=get(h,'Children');
set(a,'LineWidth',2)
xlabel('Iteration')
ylabel('Convergence Quantity')
title(['Convergence History, ',num2str(length(fvals)),' Function Calls, ',num2str(length(cvals(1,:))),' Constraint Calls'])
%legend('f_{high}(x_i)','||x_k-x^*||','||x_k-x_{k+1}||','Trust-Region Radius',3)
legend('f_{high}(x_i)','||c(x_k)||^2','||dL/dx||','||x_k-x^*||','Trust-Region Radius',3)

    
% Plot the final airfoils: **********************************************

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
    title(['Optimized Airfoil, L/D=',num2str(abs(f_hist(end))),' D=',num2str(c2+drag),' \alpha=',num2str(xk(1)*scale(1)),'^\circ'])
    xlabel('x/c')
    ylabel('t/c')
    a=get(h,'Children');
    set(a,'LineWidth',2,'MarkerSize',5)
end