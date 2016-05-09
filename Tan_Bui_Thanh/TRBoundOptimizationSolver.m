function [xc,fc,exitflag,histout,XC] = TRBoundOptimizationSolver...
    (xc,f,low,up,option)
% Tan Bui, Jan 9, 2006
% The Trust Region Interior Reflective (with Newton, QuasiNewton)
% and the Line Search Algorithm

% The Trust Region Interior Reflective Newton-CG
if strcmpi(option.OptMethod,'TIRNCG')
    [xc,fc,exitflag,histout,XC] = TIRNCG(xc,f,low,up,option);
elseif strcmpi(option.OptMethod,'TRNCG')
    [xc,fc,exitflag,histout,XC] = TRNCG(xc,f,option);
elseif strcmpi(option.OptMethod,'TIRNCGShapeCos')
    [xc,fc,exitflag,histout,XC] = TIRNCGShapeCos(xc,f,low,up,option);
    % The subspace Trust Region Interior Reflective Newton-CG
elseif strcmpi(option.OptMethod,'STIRNCG')
    [xc,fc,exitflag,histout,XC] = STIRNCG(xc,f,low,up,option);
    % The subspace Trust Region Interior Reflective Newton-CG
    % especially for 2D shape inverse EM scattering problem
elseif strcmpi(option.OptMethod,'STIRNCGShape')
    [xc,fc,exitflag,histout,XC] = STIRNCGShape(xc,f,low,up,option);
    % especially for 2D shape inverse EM scattering problem with cosines
elseif strcmpi(option.OptMethod,'STIRNCGShapeCos')
    [xc,fc,exitflag,histout,XC] = STIRNCGShapeCos(xc,f,low,up,option);
     % unscaled trust-region problem
elseif strcmpi(option.OptMethod,'TIRNCGunscalded')
    [xc,fc,exitflag,histout,XC] = TIRNCGunscalded(xc,f,low,up,option);
% The Trust Region Interior Reflective Quasi-Newton-CG
elseif strcmpi(option.OptMethod,'TIRQNCG')
    [xc,fc,exitflag,histout,XC] = TIRQNCG(xc,f,low,up,option);
% The Line Search method with Newton-CG
elseif strcmpi(option.OptMethod,'PASIPNCG')
    [xc,fc,exitflag,histout,XC] = PASIPNCG(xc,f,low,up,option);
% The Line Search method with Quasi-Newton-CG
elseif strcmpi(option.OptMethod,'PASIPQNCG')
    [xc,fc,exitflag,histout,XC] = PASIPQNCG(xc,f,low,up,option);
else, error('unsupported Optimization method, please add'); end
