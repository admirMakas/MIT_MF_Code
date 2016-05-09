% function [x,top,bot]=GenFoilFour(DV,makeplots)
%
%   Generates an airfoil using Fourier Series coefficients for the
%   camberline and the thickness. Airfoil will have sharp leading and
%   trailing edges as required for the panel method and shock-expansion
%   theory.
%
%       DV  - can be any 2n+1 length vector, and GenFoilFour will interpret
%       the vector as angle of attack, n camberline coefficients,
%       n thickness coefficients
%       makeplots - boolean, whether or not to plot the airfoil
%
%
%       x - x/c points on the airfoil
%       top/bot -  the corresponding y/c coordinate of the upper and lower
%       airfoil surfaces
%
function [x,top,bot]=GenFoilFour(DV,makeplots)
nCterms=floor(length(DV)/2);
nTterms=nCterms;
npts=50;
% camb_max=0.05;
% thick_max=0.05;
% Use either...
CambWeights=DV(2:1+nCterms);
ThickWeights=DV(2+nCterms:nCterms+nTterms+1);
% or  these...
% rand('twister',sum(100*clock))
% CambWeights=-camb_max+2*camb_max*rand(nCterms,1);
% ThickWeights=thick_max*rand(nTterms,1);
% done
x=linspace(0,1,npts);
camb=zeros(size(x));
thick=zeros(size(x));
for i=1:nCterms
    camb=camb+CambWeights(i)*sin(i*pi*x);
end
% iter=0;
% while(sum(ThickWeights(2:end))>ThickWeights(1) && iter<10)
%     ThickWeights(2:end)=ThickWeights(2:end)/2;
%     iter=iter+1;
% end

for i=1:nTterms
    thick=thick+ThickWeights(i)*sin(i*pi*x);
end
thick(thick<0)=0;
top=camb+thick;
bot=camb-thick;
if(any(isnan(top)) || any(isnan(bot)))
    makeplots=true;
    disp('problem')
end
if(makeplots)
    plot(x,camb,'g')
    hold on;
    plot(x,top,'r')
    plot(x,bot,'r')
    hold off;
    xlabel('x')
    ylabel('Airfoil Profile')
    title('Airfoil Shape')
end