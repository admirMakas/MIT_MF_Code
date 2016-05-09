% function [x,top,bot]=GenFoilPts(DV,makeplots)
%
%   Generates an airfoil using points on the upper and lower surface of the
%   airfoil. If more points are requested than points included than in the
%   design vector, the points will be linearly interpolated. The airfoil
%   will have sharp leading and trailing edges as required for the panel
%   method and shock-expansion theory.
%
%       DV  - can be any 2n+1 length vector, and GenFoilPts will interpret
%       the vector as angle of attack, n upper surface points,
%       n lower surface points (t/c)
%       makeplots - boolean, whether or not to plot the airfoil
%
%
%       x - x/c points on the airfoil
%       top/bot -  the corresponding y/c coordinate of the upper and lower
%       airfoil surfaces
%
function [x,top,bot]=GenFoilPts(DV,makeplots)
% assignin('base','DVGen',DV)
npts=250;
ntop=floor(length(DV)/2);
nbot=ntop;
x=linspace(0,1,ntop+2)';
% Use either...
top=[0,DV(2:1+ntop),0];
bot=[0,DV(2+ntop:ntop+nbot+1),0];
if(floor(npts/2)>ntop+2)
    x2=x;
    x=linspace(0,1,floor(npts/2));
    top=interp1(x2,top,x,'linear','extrap');
    bot=interp1(x2,bot,x,'linear','extrap');
end
if(any(isnan(top)) || any(isnan(bot)))
    makeplots=true;
    disp('problem')
end
if(makeplots)
    plot(x,top,'r')
    hold on;
    plot(x,bot,'m')
    hold off;
    xlabel('x')
    ylabel('Airfoil Profile')
    title('Airfoil Shape')
end