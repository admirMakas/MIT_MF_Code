close all;
clc
font=12;

fn=@f;
constr=@multiFconst;
LB=[-5;-5];
UB=[5;5];
n=length(LB);
ndiv=0;
count=1;

fid=[-1,2];
fn_param={};
constr_param={};

Xo=LB+(UB-LB).*rand(n,1);
LHS=lhsdesign(n,ndiv);
points=zeros(n,ndiv);
c1=zeros(ndiv+1,1);
c2=zeros(ndiv+1,1);
for i=1:ndiv+1
    if(i==1)
        points(:,i)=Xo;
    else
        points(:,i)=LB+(UB-LB).*LHS(:,i-1);
    end
    c1(i)=constr(points(:,i),fid(1),constr_param{:});
    c2(i)=constr(points(:,i),fid(2),constr_param{:});  
end
cvals=c2-c1;
xk=Xo;
fk=f(Xo);

pmax=50;
eta0=0;
eta1=.5;
eta2=1.25;
gamma0=1/2;
gamma1=2;
delta0=max(10,norm(xk,inf));
deltamax=10^3*delta0;
epsilon=1*10^(-4);
kd=10^-4;
alpha=.9;
mu=100;
beta=3;
theta1=10^-3;
theta2=10^-4;
theta3=10;
theta4=max(sqrt(n),10);
iter=0;
epsOpt=epsilon/100;
epsCon=epsilon/100;
options0=optimset('MaxIter',100,'Algorithm','active-set','Display','off','TolFun',epsOpt,'TolX',epsOpt,'TolCon',epsCon);
options=optimset('MaxIter',100,'Algorithm','active-set','Display','off','TolFun',epsOpt,'TolX',epsOpt,'TolCon',epsCon);
options2=optimset('MaxIter',100,'Algorithm','interior-point','Display','off','TolFun',epsOpt,'TolX',epsOpt,'TolCon',epsCon);
delta=delta0;
fd=10^-6; %epsilon/10;
rpen=exp(iter/10);

xk=Xo;
x_hist=xk;
f_hist=fk;
d_hist=delta;


corr='Gauss';
gammaR=2;
betaR=[];
params=[betaR,gammaR];
% corr='Cubic';
% params=3;
[Y2,Coeffs,linear,points,cvals,params]=RBFModel1(constr,xk,c2-c1,points,cvals,theta1,theta2,theta3,theta4,delta,deltamax,pmax,corr,params,constr_param,fid);

s=ones(n,1);

while(~(delta<=epsilon || (0) ) && iter<100)
    
% CRITICALITY CHECK!!!!
    
    % find direction sk
    s0=zeros(n,1);
    [s,fks,flag,output,lambda,grad]=fmincon(@(s) fn(xk+s,fn_param{:}),s0,[],[],[],[],LB,UB,@(s) RBFconst(xk+s,s,delta,constr,Y2,Coeffs,corr,params,constr_param,fid(1)),options0);
%     if(any(isnan(s)))
%         'interior point'
%         [s,mks,exitflag,output,lambda1,dm_dx]=fmincon(@(s) evalRBF(xk+s,fn,Y2,Coeffs,corr,params,fn_param,fid),s0,[],[],[],[],[],[],@(s) mycons(xk+s,xk,delta,constr,fn_param),options2);
%     end
    lambda1=lambda.ineqnonlin;
    if(any(isnan(s)))
        s=zeros(n,1);
    end
    
    if(iter==0)
        s_hist=s;
        lambda_hist=lambda1(1);
    else
        s_hist=[s_hist,s];
        lambda_hist=[lambda_hist,lambda1(1)];
    end
     
    
%*************************************************************
    
    % update trust region radius
    if(lambda1(1)>epsCon)
%         delta=min(gamma1*delta,deltamax);
        delta=gamma0*delta;
    elseif(lambda1(2)>epsCon) 
        delta=min(gamma1*delta,deltamax);
%         delta=gamma0*delta;
    else
        delta=gamma0*delta;
    end
    
    % update the iterate xk, based on if still feasible after update
      % new mks value
    i=1;
    found=-1;
    while(i<=length(points(1,:)) && found==-1)
        if(norm(points(:,i)-(xk+s))<=1e-12)
            found=i;
        end
        i=i+1;
    end
    if(found==-1)
        c1s=constr(xk+s,fid(1),constr_param{:});
        c2s=constr(xk+s,fid(2),constr_param{:});
        points=[points,xk+s];
        cvals=[cvals;c2s-c1s];
    else
        c1s=constr(xk+s,fid(1),constr_param{:});
        c2s=c1s+cvals(found);
    end
    
%     [Y2,Coeffs,linear,points,cvals,params]=RBFModel1(constr,xk+s,c2s-c1s,points,cvals,theta1,theta2,theta3,theta4,delta,deltamax,pmax,corr,params,constr_param,fid);
%     c2s_est=evalRBF(xk+s,constr,Y2,Coeffs,corr,params,constr_param,fid);
%     if(abs(c2s_est-c2s)>epsCon)
%         c2s_est
%         c2s
%     end
    
    xk_old=xk;
    if(c2s<=0 || c2s<=c2)
        xk=xk+s;
        fk=fks;
        c2=c2s;
        c1=c1s;
    else
        xk=xk;
        fk=fk;
    end
    
    
    [Y2,Coeffs,linear,points,cvals,params]=RBFModel1(constr,xk,c2-c1,points,cvals,theta1,theta2,theta3,theta4,delta,deltamax,pmax,corr,params,constr_param,fid);
    c2_est=evalRBF(xk,constr,Y2,Coeffs,corr,params,constr_param,fid);
    if(abs(c2_est-c2)>epsCon)
        c2s_est
        c2
    end
    iter=iter+1;
    
%     disp(['iter: ',num2str(iter)]);
    x_hist=[x_hist,xk];
    f_hist=[f_hist,fk];
    d_hist=[d_hist,delta];
    
%     if(delta<epsilon/1000)
%         delta=epsilon/1000;
%     end
%     rpen=exp(iter/10);
end
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
subplot(1,2,1)
h=gca;
set(h,'FontSize',font)
semilogy(abs(f_hist))
hold on;
semilogy(x_hist2,'k')
semilogy(dx_hist2,'r')
semilogy(abs(d_hist),'g')

hold off;
xlabel('Iteration')
ylabel('Convergence Quantity')
title('Convergence History')
legend('f_{high}(x_i)','||x_k-x^*||','||x_k-x_{k+1}||','Trust-Region Radius',3)
subplot(1,2,2)
h=gca;
set(h,'FontSize',font)
plot(x_hist(1,:),x_hist(2,:))
hold on;
plot(x_hist(1,1),x_hist(2,1),'ob')
plot(x_hist(1,end),x_hist(2,end),'*r')
% thetas=linspace(0,2*pi,50);
% plot(2*cos(thetas),2*sin(thetas),':k')
val=1000;
text(x_hist(1,end),x_hist(2,end),['(',num2str(round(val*x_hist(1,end))/val),',',...
    num2str(round(val*x_hist(2,end))/val),')'])
hold off;
xlabel('x_1')
ylabel('x_2')
title('Optimization Path')

figure
h=gca;
set(h,'FontSize',font)
semilogy(abs(f_hist))
hold on;
semilogy(x_hist2,'k')
semilogy(dx_hist2,'r')
semilogy(abs(d_hist),'g')
% semilogy(s2,':k')
% semilogy(dx_hist2,'c')
semilogy(abs(lambda_hist-lambda_hist(end)),'m')
hold off;
a=get(h,'Children');
set(a,'LineWidth',2)
xlabel('Iteration')
ylabel('Convergence Quantity')
title('Convergence History')
legend('f_{high}(x_i)','||x_k-x^*||','||x_k-x_{k+1}||','Trust-Region Radius',3)
% legend('f_{high}(x_i)','|f_{high}(x_i)-f_{high}(x^*)|','Trust-Region
% Radius',3)


% options2=optimset('Algorithm','interior-point','Display','off','TolFun',epsOpt,'TolX',epsOpt,'TolCon',epsOpt);
count=1;
[x2,ftrue,exitflag,output,lambda_true,df_dx]=fmincon(@(x) fn(x) ,points(:,1),[],[],[],[],LB,UB,@(x) constr(x,fid(2),constr_param{:}),options2);
% output.funcCount
xk
SQP_ans=x2
disp(['Number of high-fidelity calls=',num2str(length(points)),'.'])
disp(['Number of SQP constraint calls=',num2str(count)'.'])
if(norm(x2-xk)>epsilon)
    disp('Bad solution, or different solution from fmincon direct solution')
    norm(xk)
end

% options=optimset('MaxFunEvals',3000,'Algorithm','active-set','Display','iter','TolFun',epsOpt,'TolX',epsOpt,'TolCon',epsOpt);
% [x_best,mks,exitflag,output,lambda1,dm_dx]=fmincon(@(xk)
% analyze(xk,2,M,gamma,FoilGen),Xo,[],[],[],[],[],[],@(xk) myAFconst(xk,fn_param{:}),options);

% step=max(delta,1e-5);
% npts=20;
% xvals=linspace(xk(1)-step,xk(1)+step,npts);
% yvals=linspace(xk(2)-step,xk(2)+step,npts);
% [X1,Y1]=meshgrid(xvals,yvals);
% Z1=zeros(npts,npts);
% for i=1:npts
%     for j=1:npts
% %         Z1(i,j)=fn([X1(i,j);Y1(i,j)],0);
%         Z1(i,j)=evalRBF([X1(i,j);Y1(i,j)],fn,Y2,Coeffs,corr,params);
%     end
% end
% index=[];
% for i=1:length(points(1,:))
%     if(points(1,i)>=xk(1)-step && points(1,i)<=xk(1)+step && ...
%             points(2,i)>=xk(2)-step && points(2,i)<=xk(2)+step)
%         index=[index,i];
%     end
% end
% figure
% subplot(1,2,1)
% surf(X1,Y1,Z1)
% subplot(1,2,2)
% % plot([xk(1)-step,xk(1)-step,xk(1)+step,xk(1)+step,xk(1)-step],...
% %     [xk(2)-step,xk(2)+step,xk(2)+step,xk(2)-step,xk(2)-step],':k')
% thetas=linspace(0,2*pi,25);
% plot(xk(1)+step*cos(thetas),xk(2)+step*sin(thetas),':k')
% hold on;
% plot(points(1,index),points(2,index),'or')
% hold off;
% axis equal
% 
% npts=20;
% xvals=linspace(LB(1),UB(1),npts);
% yvals=linspace(LB(2),UB(2),npts);
% [X1,Y1]=meshgrid(xvals,yvals);
% Z1=zeros(npts,npts);
% Z2=zeros(npts,npts);
% for i=1:npts
%     for j=1:npts
%         Z1(i,j)=fn([X1(i,j);Y1(i,j)],-1);
%         Z2(i,j)=fn([X1(i,j);Y1(i,j)],2);
%     end
% end
% figure
% subplot(1,2,1)
% h=gca;
% set(h,'FontSize',font)
% surf(X1,Y1,Z1)
% xlabel('x_1')
% ylabel('x_2')
% title('Low-Fidelity Model')
% subplot(1,2,2)
% h=gca;
% set(h,'FontSize',font)
% surf(X1,Y1,Z2)
% xlabel('x_1')
% ylabel('x_2')
% title('High-Fidelity Model')