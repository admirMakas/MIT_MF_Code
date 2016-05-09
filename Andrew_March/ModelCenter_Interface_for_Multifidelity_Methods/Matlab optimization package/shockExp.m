% function [Cl,Cd]=shockExp(M,alpha,gamma,x,top,bot)
%
%   Solves for the supersonic lift and drag coefficents for an airfoil with
%   a sharp leading and trailing edge using shock-expansion theory.
%   Shock-expansion theory will fail if there is a strong shockwave
%   (subsonic downstream Mach number), or if the shockwave is detached. In
%   either case, when detected the resulting output will be Cl=0, Cd=-1.
%
%       M_inf - freestream Mach number
%       alpha - freestream angle of attack
%       gamma - ratio of specific heats
%       x - vector of x/c points
%       top - vector of points on the airfoil upper surface y/c
%       bot - vector of points on the airfoil lower surface y/c
%
%       Cl - lift coefficient
%       Cd - drag coefficient
%
function [Cl,Cd]=shockExp(M_inf,alpha,gamma,x,top,bot)
% Change plot ops to 
%   (1) plot the cp distribution? 
%   (2) in a new figure-true, subplot for false
plotops=[0,0];
angles=zeros(2,length(top)-1);
alpha=alpha*pi/180;
x2=zeros(1,length(x)-1);
lens=angles;
Cp=angles;
P=zeros(2,length(top));
P(:,1)=[1;1];
d_angles=angles;
Cl=0;
Cd=0;
M_top=M_inf;
M_bot=M_inf;
i=1;
failed=false;
eps=1e-14;
while i<=length(x2) && ~failed
    angles(1,i)=atan2(top(i+1)-top(i),x(i+1)-x(i))-alpha;
    angles(2,i)=-atan2(bot(i+1)-bot(i),x(i+1)-x(i))+alpha;
    lens(1,i)=sqrt((x(i+1)-x(i))^2+(top(i+1)-top(i))^2);
    lens(2,i)=sqrt((x(i+1)-x(i))^2+(bot(i+1)-bot(i))^2);
    x2(i)=(x(i+1)+x(i))/2;
    if(i>1)
        d_angles(:,i)=angles(:,i)-angles(:,i-1);
    else
        d_angles(1,i)=angles(1,i);
        d_angles(2,i)=angles(2,i);
    end
    if(any(abs(d_angles(:,i))>=pi/3))
        d_angles(:,i)=zeros(2,1);
        failed=true;
    end
    if(abs(d_angles(1,i))<=eps && i>1)
        % no angle change;
        P(1,i+1)=P(1,i);
    elseif(d_angles(1,i)>0)
        % Compress top
        M=M_top;
        theta=d_angles(1,i);
        a=1;
        b=-(M^2+2)/M^2-gamma*sin(theta)^2;
        c=(2*M^2+1)/M^4+((gamma+1)^2/4+(gamma-1)/M^2)*sin(theta)^2;
        d=-cos(theta)^2/M^4;
        vals=roots([a,b,c,d]);
        betas=asin(sqrt(vals));
        beta=0;
        if(isreal(betas))
            betas=sort(betas);
            beta=betas(2);
%             beta2=fzero(@(beta) Oblique(beta,theta,M,gamma),asin(1/M));
        else
            failed=true;
        end
        Mn1=M*sin(beta);
        M_top=sqrt((1+((gamma-1)/2)*Mn1^2)/(gamma*Mn1^2-(gamma-1)/2))/(sin(beta-theta));
        P(1,i+1)=(1+2*gamma/(gamma+1)*(Mn1^2-1))*P(1,i);
        if(beta<asin(1/M) || M_top<1 || beta>89*pi/180)
            failed=true;
%             betas
%             beta
%             M
%             theta
        end
    else
        % expand top
        M=M_top;
        v=Prandtl(M,gamma); % find P-M angle
        v=v-d_angles(1,i); % Find angle of expanded flow
%         M_top=fzero(@(M) Prandtl(M,gamma)-v,M*1.1);
        M_top=fzero(@(M) Prandtl(M,gamma)-v,[M,5*M]);
        P(1,i+1)=(((1+(gamma-1)/2*M^2)/(1+(gamma-1)/2*M_top^2))^(gamma/(gamma-1)))*P(1,i);
    end
    if(abs(d_angles(2,i))<=eps && i>1)
        % no angle change;
        P(2,i+1)=P(2,i);
    elseif(d_angles(2,i)>0)
        % Compress bot
        M=M_bot;
        theta=d_angles(2,i);
        a=1;
        b=-(M^2+2)/M^2-gamma*sin(theta)^2;
        c=(2*M^2+1)/M^4+((gamma+1)^2/4+(gamma-1)/M^2)*sin(theta)^2;
        d=-cos(theta)^2/M^4;
        vals=roots([a,b,c,d]);
        betas=asin(sqrt(vals));
        beta=0;
        if(isreal(betas))
            betas=sort(betas);
            beta=betas(2);
%             beta2=fzero(@(beta) Oblique(beta,theta,M,gamma),asin(1/M));
%             if(abs(beta-beta2)>=10^(-4))
%                 beta
%                 betas
%                 M
%                 theta
%             end
        else
            failed=true;
        end
        Mn1=M*sin(beta);
        M_bot=sqrt((1+((gamma-1)/2)*Mn1^2)/(gamma*Mn1^2-(gamma-1)/2))/(sin(beta-theta));
        P(2,i+1)=(1+2*gamma/(gamma+1)*(Mn1^2-1))*P(2,i);
        
        if(beta<asin(1/M) || M_bot<1 || beta>89*pi/180)
            failed=true;
%             betas
%             beta
%             M
%             theta
        end
        
    else
        % expand bot
        M=M_bot;
        v=Prandtl(M,gamma); % find P-M angle
        v=v-d_angles(2,i); % Find angle of expanded flow
%         M_bot=fzero(@(M) Prandtl(M,gamma)-v,M*1.1);
        M_bot=fzero(@(M) Prandtl(M,gamma)-v,[M,5*M]);
        P(2,i+1)=(((1+(gamma-1)/2*M^2)/(1+(gamma-1)/2*M_bot^2))^(gamma/(gamma-1)))*P(2,i);
    end
    Cp(1,i)=(P(1,i+1)-1)/(gamma/2*M_inf^2);
    Cp(2,i)=(P(2,i+1)-1)/(gamma/2*M_inf^2);
    Cl=Cl+Cp(2,i)*cos(angles(2,i))*lens(2,i)-...
        Cp(1,i)*cos(angles(1,i))*lens(1,i);
    Cd=Cd+Cp(2,i)*sin(angles(2,i))*lens(2,i)+...
        Cp(1,i)*sin(angles(1,i))*lens(1,i);
    i=i+1;
    if(M_top<1 || M_top<1)
        failed=true;
    end
end
if(plotops(1))
    if(plotops(2))
        figure
    else
        subplot(1,3,3)
    end
    plot(x2,Cp)
    hold on;
    plot(x,[top;bot],'r')
    hold off;
    legend('Cp-top','Cp-bot')
    title('Shock Exp')
end
if(failed)
    Cd=-1;
%     Cd=10+max(max(d_angles(1,:)),max(d_angles(2,:)));
    Cl=0;
%     M_top
%     M_bot
end