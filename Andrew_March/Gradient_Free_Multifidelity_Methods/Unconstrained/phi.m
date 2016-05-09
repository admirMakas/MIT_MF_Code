% function phi_ij=phi(r,corr,varargin)
%
% Computes the radial basis correlation for all the twice continuously
% differentiable radial basis functions in Wild, Shoemaker "Global
% Convergence of Radial Basis Function Trust-Region Algorithms"
%                       phi(r)              Order   Parameters
%   Cubic:              r^beta              2       beta (2,4)
%   Multiquadratic I:   (gamma^2+r^2)^beta  2       beta (1,2), gamma>0
%   Multiquadratic II:  -(gamma^2+r^2)^beta 1       beta (0,1), gamma>0
%   Inv. Multiquadratic: (gamma^2+r^2)^-beta 0      beta >0, gamma>0
%   Gaussian:           exp^(-r^2/gamma^2)  0       gamma>0

function phi_ij=phi(r,corr,params)
switch corr
    case {'cubic','Cubic'}
        if(params(1)<=2 || params(1)>=4)
            warndlg('Invalid Parameters for Cubic Correlation.','phi');
        end
        phi_ij=r^params(1);
    case {'MultiI','multiI','Multiquadratic I','multiquadratic I','MultiquadraticI','multiquadraticI'}
        phi_ij=(params(2)^2+r^2)^params(1);
        if(params(1)<=1 || params(1)>=2 || params(2)<=0)
            warndlg('Invalid Parameters for Multiquadratic I Correlation.','phi');
        end
    case {'MultiII','multiII','Multiquadratic II','multiquadratic II','MultiquadraticII','multiquadraticII'}
        phi_ij=-(params(2)^2+r^2)^params(1);
        if(params(1)<=0 || params(1)>=1 || params(2)<=0)
            warndlg('Invalid Parameters for Multiquadratic I Correlation.','phi');
        end  
    case {'InvMulti','invmulti','Inv. Multiquadratic','inv. multiquadratic','InvMultiquadratic','invmultiquadratic'}
        phi_ij=(params(2)^2+r^2)^(-params(1));
        if(params(1)<=0 || params(2)<=0)
            warndlg('Invalid Parameters for Multiquadratic I Correlation.','phi');
        end 
    case {'Gaussian','Gauss','gaussian','gauss'}
        if(params(1)<=0)
            warndlg('Invalid Parameters for Cubic Correlation.','phi');
        end
        phi_ij=exp(-r^2/params(1)^2);
end