function [F,J] = SubTRCost(z,varargin)
% Tan Bui: Jan 7, 2007
% Cost functional for the Lagrangian
    
    load SubTRdata
    
    del = (z*c(2))^2-c(3)*(z^2*c(1)-Delta^2);
    alpha1 = ((-z*c(2))+sqrt(del))/c(3);
    alpha2 = ((-z*c(2))-sqrt(del))/c(3);
    
% $$$     if ~(isreal(alpha1))||~(isreal(alpha2)),
% $$$         alpha1, alpha2
% $$$         error('alpha must be real')
% $$$     end
    
    z(2,1) = alpha1; 
    m1 = z(1)*b(1)+z(2)*b(2)+0.5*z(1)^2*a(1)+z(1)*z(2)*a(2)+0.5* ...
        z(2)^2*a(3);
    z(2,1) = alpha2;
    m2 = z(1)*b(1)+z(2)*b(2)+0.5*z(1)^2*a(1)+z(1)*z(2)*a(2)+0.5* ...
        z(2)^2*a(3);
    
    if (m1 > m2)
        z(2,1) = alpha2;
        F = m2;
    else
        z(2,1) = alpha1;
        F = m1;
    end
    
    save SubTRSolution z
    
    if nargout > 1,
        z(3) = -(b(2)+z(1)*a(2)+z(2)*a(3))/(2*z(1)*c(2)+2*z(2)* ...
                                             c(3));
        J =  b(1)+z(1)*a(1)+z(2)*a(2)+2*z(3)*z(1)*c(1)+...
              2*z(3)*z(2)*c(2);
        
        save SubTRHessdata z
    end
    
    
    
    
    
    
    
    
    
    
    
    
% $$$     m = z(1)*b(1)+z(2)*b(2)+0.5*z(1)^2*a(1)+z(1)*z(2)*a(2)+0.5* ...
% $$$         z(2)^2*a(3);
% $$$     
% $$$     con = z(1)^2*c(1)+2*z(1)*z(2)*c(2)+z(2)^2*c(3)- ...
% $$$                   Delta^2;
% $$$     F = m + z(3)*con;
% $$$     
% $$$     if nargout > 1
% $$$        J = zeros(3,1);
% $$$        J(1) = b(1)+z(1)*a(1)+z(2)*a(2)+2*z(3)*z(1)*c(1)+...
% $$$               2*z(3)*z(2)*c(2);
% $$$        J(2) = b(2)+z(1)*a(2)+z(2)*a(3)+2*z(3)*z(1)*c(2)+...
% $$$               2*z(3)*z(2)*c(3);
% $$$        J(3) = con;
% $$$     end