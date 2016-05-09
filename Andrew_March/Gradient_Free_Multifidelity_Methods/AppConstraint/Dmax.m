% function [D]=Dmax(DV,fid,M,gamma,FoilGen,ncpts,scale)
% 
%   This is the function that chooses the proper analysis to run based on
%   the analysis fidelity level and airfoil parametrization used.
%
%       Within the function, the variable L_D_or_D needs to be set. If the
%       value is 1 the function will return the L/D ratio; if it is zero
%       the function will return the drag coefficient.
%       DV - the design vector, alpha+2N terms
%       fid - the fidelity level:
%           1 - Thin airfoil theory
%           2 - linear supersonic panel method
%           3 - shock-expansion theory
%           4 - Cart3D
%       M - Mach number
%       gamma - specific heat ratio
%       FoilGen - parameter indicating the parameterization used
%           0 - Fourier series coefficients
%           1 - Spline points
%           2 - Points on the airfoil
%       ncpts - number of control spline points on the airfoil
%       scale - scale vector to scale the design variables
%
%       D - if L_D_or_D=0, the drag coefficient, if L_D_or_D=1 the lift to
%       drag ratio
%
function D=Dmax(DV,fid,drag,M,gamma,FoilGen,ncpts,scale)
% evalin('base','count=count+1;')
% DV(1)=DV(1)*10;
[dim1,dim2]=size(DV);
if(dim2>dim1)
    DV=DV.*scale';
else
    DV=DV.*scale;
end
L_D_or_D=0;
if(L_D_or_D)
    % FoilGen=1; % 0 for Fourier Series, 1 for Spline, 2 for points on surface
    switch FoilGen
        case 0
            [x,top,bot]=GenFoilFour(DV,false);
        case 1
            [x,top,bot]=GenFoilSpline(DV,false);
        case 2
            [x,top,bot]=GenFoilPts(DV,false);
    end
    top=top';
    bot=bot';
    alpha=DV(1);
    if(any(isnan(DV)))
        D=1;
    else
        switch fid
            case 1
                camb_surf=mean([top;bot],1);
                [CL,CD]=thinAirfoil(M,alpha,x,camb_surf);
            case 2
                [CL,CD]=thickAirfoil(M,alpha,x,top,bot);
            case 3
                [CL,CD]=shockExp(M,alpha,gamma,x,top,bot);
                if(CD<0)
                    [CL,CD]=thickAirfoil(M,alpha,x,top,bot);
                end
            case 4
                [CL,CD]=CartMesh(M,alpha,x,top,bot,'Components.i');
                if(CD<0)
                    [CL,CD]=thickAirfoil(M,alpha,x,top,bot);
                end

        end
        if(CD<=0)
            D=1;
        else
            D=CD;
        end
        D=-CL/D;
    end
else
    % function D=analyze(DV,fid,M,gamma,FoilGen)
    % FoilGen=1; % 0 for Fourier Series, 1 for Spline, 2 for points on surface
    switch FoilGen
        case 0
            [x,top,bot]=GenFoilFour(DV,false);
        case 1
            [x,top,bot]=GenFoilSpline(DV,false);
        case 2
            [x,top,bot]=GenFoilPts(DV,false);
    end
    top=top';
    bot=bot';
    alpha=DV(1);
    if(any(isnan(DV)))
        D=100;
    else
        switch fid
            case 1
                camb_surf=mean([top;bot],1);
                [CL,CD]=thinAirfoil(M,alpha,x,camb_surf);
            case 2
                [CL,CD]=thickAirfoil(M,alpha,x,top,bot);
    %             [CL,CD]=thickAirfoil2(M,alpha,gamma,x,top,bot);
            case 3
                [CL,CD]=shockExp(M,alpha,gamma,x,top,bot);
                if(CD<0)
                    [CL,CD]=thickAirfoil(M,alpha,x,top,bot);
                end
            case 4
                [CL,CD]=CartMesh(M,alpha,x,top,bot,'Components.i');
                if(CD<0)
                    [CL,CD]=thickAirfoil(M,alpha,x,top,bot);
                end

        end
        if(CD<=0)
            D=100;
        else
            D=CD;
        end
%         thick=top-bot;
%         if(0.05-max(thick)>0) % garentee thickness
%             D=D+1000*(0.05-max(thick))^2;
%         end
%         if(-min(thick(2:length(thick)-1))>0) % all positive thickness
%             D=D+10*(min(thick(2:length(thick)-1)))^2;
%         end
% %         if(-pi/2+(atan2(top(2)-top(1),x(2)-x(1))-DV(1)*pi/180)>0)
% %             D=D+10*(-pi/2+(atan2(top(2)-top(1),x(2)-x(1))-DV(1)*pi/180))^2;
% %         end 
        
    end
    D=(D-drag);
end