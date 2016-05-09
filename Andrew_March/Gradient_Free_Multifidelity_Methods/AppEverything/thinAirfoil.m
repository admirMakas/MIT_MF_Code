% function [Cl,Cd]=thinAirfoil(M,alpha,x,camb)
%
%   Collapses an airfoil to its camberline and solves for the supersonic
%   lift and drag coefficents assuming linearized small-disturbance  flow.
%
%       M - freestream Mach number
%       alpha - freestream angle of attack
%       x - vector of x/c points
%       camb - vector of points on the camberline y/c
%
%       Cl - lift coefficient
%       Cd - drag coefficient
%
function [Cl,Cd]=thinAirfoil(M,alpha,x,camb)
plotops=[0,0];
angles=zeros(1,length(camb)-1);
alpha=alpha*pi/180;
x2=angles;
Cp=angles;
lens=angles;
d_angles=angles;
Cl=0;
Cd=0;
for i=1:length(camb)-1
    angles(i)=atan2(camb(i+1)-camb(i),x(i+1)-x(i))+alpha;
    lens(i)=sqrt((x(i+1)-x(i))^2+(camb(i+1)-camb(i))^2);
    x2(i)=(x(i+1)+x(i))/2;
    if(i>1)
        d_angles(i)=angles(i)-angles(i-1);
        Cp(i)=Cp(i-1)+2*d_angles(i)/sqrt(M^2-1);
    else
        d_angles(i)=angles(i);
        Cp(i)=2*d_angles(i)/sqrt(M^2-1);
    end
    Cl=Cl+2*Cp(i)*cos(angles(i))*lens(i);
    Cd=Cd+2*Cp(i)*sin(angles(i))*lens(i);
end
if(plotops(1))
    if(plotops(2))
        figure
    else
        subplot(1,3,1)
    end
    plot(x2,[-1*Cp;Cp])
    hold on;
    plot(x,camb,'r')
    hold off;
    legend('Cp-top','Cp-bot')
    title('Thin Approx')
end