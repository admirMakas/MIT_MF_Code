% function [x,top,bot]=GenFoilSpline(DV,makeplots)
%
%   Generates an airfoil using spline points on the upper and lower surface
%   of the airfoil. If more points are requested than points included than
%   in the design vector, the points will be interpolated using a cubic
%   spline. The airfoil will have sharp leading and trailing edges as
%   required for the panel method and shock-expansion theory.
%
%       DV  - can be any 2n+1 length vector, and GenFoilSpline will interpret
%       the vector as angle of attack, n upper surface points,
%       n lower surface points (t/c)
%       makeplots - boolean, whether or not to plot the airfoil
%
%       x - x/c points on the airfoil
%       top/bot -  the corresponding y/c coordinate of the upper and lower
%       airfoil surfaces
%
% NOTE: You can change the npts varable which is the number of surface
% points.
function [x,top,bot]=GenFoilSpline(DV,makeplots)
[m,n]=size(DV);
if(m>n)
    DV=DV';
end
% assignin('base','DVGen',DV)
ntop=floor(length(DV)/2);
nbot=ntop;
xctop=linspace(0,1,ntop+2)';
xcbot=linspace(0,1,nbot+2)';
npts=350; % Feel free to change....
% Use either...
Tpts=DV(2:1+ntop);
Bpts=DV(2+ntop:ntop+nbot+1);

x=linspace(0,1,npts);

topF=fit(xctop,[0;Tpts';0],'cubicspline');
botF=fit(xcbot,[0;Bpts';0],'cubicspline');
top=topF(x);
bot=botF(x);
thick=top-bot;
iter=0;
% while(any(thick<0) && iter<10)
%     top(thick<0)=top(thick<0)+1e-5;
%     thick=top-bot;
%     iter=iter+1;
% end
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