% function [Cl,Cd]=thickAirfoil(M_inf,alpha,x,top,bot)
%
%   Solves for the supersonic lift and drag coefficents for an arbitrary
%   airfoil using linearized small-disturbance theory.
%
%       M_inf - freestream Mach number
%       alpha - freestream angle of attack
%       x - vector of x/c points
%       top - vector of points on the airfoil upper surface y/c
%       bot - vector of points on the airfoil lower surface y/c
%
%       Cl - lift coefficient
%       Cd - drag coefficient
%
function [Cl,Cd]=thickAirfoil(M_inf,alpha,x,top,bot)
% Change plot ops to 
%   (1) plot the cp distribution? 
%   (2) in a new figure-true, subplot for false
plotops=[0,0];

angles=zeros(2,length(top)-1);
alpha=alpha*pi/180;
x2=zeros(1,length(x)-1);
lens=angles;
Cp=angles;
d_angles=angles;
Cl=0;
Cd=0;
M_top=M_inf;
M_bot=M_inf;
for i=1:length(x2)
    angles(1,i)=atan2(top(i+1)-top(i),x(i+1)-x(i))-alpha;
    angles(2,i)=-atan2(bot(i+1)-bot(i),x(i+1)-x(i))+alpha;
    % need turning angle,...
%     angles(1,i)=atan((top(i+1)-top(i))/(x(i+1)-x(i)));
%     angles(2,i)=atan((bot(i+1)-bot(i))/(x(i+1)-x(i)));
    lens(1,i)=sqrt((x(i+1)-x(i))^2+(top(i+1)-top(i))^2);
    lens(2,i)=sqrt((x(i+1)-x(i))^2+(bot(i+1)-bot(i))^2);
    x2(i)=(x(i+1)+x(i))/2;
    if(i>1)
        d_angles(:,i)=angles(:,i)-angles(:,i-1);
    else
        d_angles(1,i)=angles(1,i);
        d_angles(2,i)=angles(2,i);
    end
    if(i>1)
        Cp(1,i)=Cp(1,i-1)+2*d_angles(1,i)/sqrt(M_top^2-1);
        Cp(2,i)=Cp(2,i-1)+2*d_angles(2,i)/sqrt(M_bot^2-1);
    else
        Cp(1,i)=2*d_angles(1,i)/sqrt(M_inf^2-1);
%         angles(2,i)=angles(2,i)+alpha;
        Cp(2,i)=2*d_angles(2,i)/sqrt(M_inf^2-1);
    end
%     Cl=Cl+Cp(2,i)*cos(sum(angles(2,1:i)))*lens(2,i)-...
%         Cp(1,i)*cos(sum(angles(1,1:i)))*lens(1,i);
%     Cd=Cd+Cp(2,i)*sin(sum(angles(2,1:i)))*lens(2,i)-...
%         Cp(1,i)*sin(sum(angles(1,1:i)))*lens(1,i);
    Cl=Cl+Cp(2,i)*cos(angles(2,i))*lens(2,i)-...
        Cp(1,i)*cos(angles(1,i))*lens(1,i);
    Cd=Cd+Cp(2,i)*sin(angles(2,i))*lens(2,i)+...
        Cp(1,i)*sin(angles(1,i))*lens(1,i);
    
end
if(plotops(1))
    if(plotops(2))
        figure
    else
        subplot(1,3,2)
    end
    plot(x2,Cp)
    hold on;
    plot(x,[top;bot],'r')
    hold off;
    legend('Cp-top','Cp-bot')
    title('Linear Approx')
end