% function [c,ceq]=myAFconst(DV,M,gamma,FoilGen,ncpts,scale)
%  This function is the geometric constraints for the airfoil. In this
%  version the upper and lower bounds for the airfoils are also included
%  for simplicity, although this will be inefficient with fmincon.
%
%      Constraint 1: max t/c>=5%
%      Constraint 2: min t/c>=0%
%      Constraint 3: (in general not active or important) leading edge
%      angle of attack less than 90^o.
%
%   Inputs: 
%       DV - Design vector
%       M - Mach number
%       gamma - specific heat ratio
%       FoilGen - Which airfoil generator to use
%       ncpts - number of surface control points
%       scale - scale vector used to scale design variables
%
%   Outputs:
%       c - inequality constraints (fmincon format)
%       ceq - equality constraints (empty)
%
function [c,ceq]=myAFconst(DV,M,gamma,FoilGen,ncpts,scale)
% DV=DV.*scale;
% ADD LB and UB
LB=[-5;zeros(ncpts,1);zeros(ncpts,1)-.03];
UB=[5;zeros(ncpts,1)+.03;zeros(ncpts,1)];

DV=DV.*scale;
alpha=DV(1);
DV=DV';
% if smoothness constraint 4, otherwise 3
c=zeros(3,1);
switch FoilGen
    case 0
        [x,top,bot]=GenFoilFour(DV,false);
    case 1
        [x,top,bot]=GenFoilSpline(DV,false);
        top=top';
        bot=bot';
    case 2
        [x,top,bot]=GenFoilPts(DV,false);
end
thick=top-bot;
d_angles1=atan(diff(top)./diff(x));
d_angles2=atan(diff(bot)./diff(x));
d_angles1=diff(d_angles1);
d_angles2=diff(d_angles2);
% max_angle=max(abs([d_angles1,d_angles2]));

c(1)=0.05-max(thick); % garentee thickness
c(2)=-min(thick(2:length(thick)-1)); % all positive thickness
c(3)=-pi/2+(atan2(top(2)-top(1),x(2)-x(1))-alpha*pi/180);
% c(4)=max_angle-.02;

A2=eye(length(LB));
B=[-LB;UB];
c2=[-A2;A2]*(DV')-B;
c=[c;c2];
ceq=[];