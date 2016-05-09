function [xc,fc,exitflag,histout,XC] = TRNCG(xc,f,option)
% TRNCG Copyright and License information for TRNCG  
%      TRNCG is the subspace trust region interior reflective Newton-CG method
% 
%      The following information is a copy of the License file in the TRNCG
%      distribution.
%   
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   
%    TRNCG, Copyright (c) 2007 by Tan Bui
%    All Rights Reserved. Massachusetts Institute of Technology
%    
%    TRNCG License:
%    
%        Your use or distribution of TRNCG or any modified version of
%        TRNCG implies that you agree to this License.
%    
%        THIS MATERIAL IS PROVIDED AS IS, WITH ABSOLUTELY NO WARRANTY
%        EXPRESSED OR IMPLIED.  ANY USE IS AT YOUR OWN RISK.
%    
%        Permission is hereby granted to use or copy this program, provided
%        that the Copyright, this License, and the Availability of the original
%        version is retained on all copies.  User documentation of any code that
%        uses TRNCG or any modified version of TRNCG code must cite the
%        Copyright, this License, the Availability note, and "Used by permission."
%        Permission to modify the code and to distribute modified code is granted,
%        provided the Copyright, this License, and the Availability note are
%        retained, and a notice that the code was modified is included.  This
%        software was developed with support from the National Science Foundation,
%        and is provided to you free of charge.
%    
%    Availability:
%    
%        web.mit.edu/tanbui/Public/TRNCGcode/
%   
%    Main reference:
%    
%    Numerical Optimization, by Nocedal and Wright
%  
%    Bug report:
%    tanbui@alum.mit.edu or tansweet@gmail.com
% 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% xc, low and up must be column vectors
% exitflag: 1: cost below the tolerance; 2: step is less than tol;
%           3: TR size is too small; 4: Iter exceeds MaxIter
%           5: 'Scaled gradient is less than the tolerance'



% hardcoded all the parameters
mu = 1.e-1; beta = 0.25; nu = 0.75; gam0 = 0.0625; gam1 = 0.5; gam2 = 2.0;
Laml = 1; n = length(xc);

if (nargin<3),
    option.TolF  = 1.e-10; option.TolX  = sqrt(option.TolF)/10;
    option.TolGrad = option.TolX/10; option.MaxIter = 100;; 
    option.AffineScaling = 'huuIIb';
end

it = 1; histout = zeros(1,5); XC = []; 
if (nargout > 4), XC(:,end+1) = xc; end

fprintf(option.SubTRFile,'\n');
fprintf(option.SubTRFile,'Iter \t\t  F \t\t||g|| \t\t Step\t\t  TR size\t CGIter\n');

HessVect = option.HessVect;

% The main loop of the iterative optimization
while (1)

    % compute the initial function and gradient
    [fc,gc]=feval(f,xc);
    histout(end+1,1) = fc;
    % Compute the Affine scaling
    histout(end,2) = norm(gc);
 
    % initial TR radius 
    if (it == 1), Lamu = 20*norm(gc); Delta = 0.5e-2*Lamu; Deltao = Delta; end
    
    if ((it == 1)&(abs(histout(end,1)) < option.TolF)), histout(end,5)=0;
        OptimizationDisplay(histout,option); exitflag = 1; return; end
    % if scaling gradient is small, exit
    if (histout(end,2) < option.TolGrad), histout(end,3:5) = 0;
        OptimizationDisplay(histout,option); exitflag = 5; return; end

    % Inexact Newton CG 
    [s,dirs,CGIter] = NewtonCG(xc, f, gc, Delta,option.pcv,option.HessVect);
    
    pred = gc.'*s + 0.5*s'*(feval(HessVect,xc,s));
    
    if (norm(s) < option.TolX), exitflag = 2; histout(end,3) = norm(s);
        OptimizationDisplay(histout,option); return; end
        
        if (pred > 0.0), pred, error('pred cannot be positive!'); end
    
    % actual reduction versus predicted reduction
    while 1,
        fn = feval(f,xc+s); 
        if ~isreal(fn)
            s = s*0.8;
        else
            break;
        end
    end
    if (abs(fn-histout(end,1)) <= option.TolF*(1+abs(histout(end,1)))), 
            exitflag = 1; xc = xc + s; histout(end,3) = norm(s);
       histout(end,4:5) = 0; OptimizationDisplay(histout,option); return; end
        
    ared = fn-fc;
    rho = ared/pred;
   
    % If the actual reduction is reasonable, take the new point
    % and adjust the TR radius and move on
    if ((rho > mu) & (ared < 0.0)) | (ared < 0.0)
        xc = xc + s; if (nargout > 4), XC(:,end+1) = xc; end
        if (rho<= beta) & (norm(s) < gam1*Delta),
            Deltao=Delta; Delta = gam1*Delta;
        elseif rho > nu,
            if (rho > Laml)&&(norm(s) >= 0.8*Delta),
                Deltao=Delta; Delta = min(gam2*Delta,Lamu);
            else
               Deltao=Delta; Delta = min(max(Delta,gam2*norm(s)),Lamu);
            end
% $$$         elseif (rho>nu) & (norm(s) >= 0.8*Delta), 
% $$$             Deltao = Delta; Delta = min(gam2*Delta,Lamu);
        else, Deltao = Delta; end

    else
        fprintf(option.SubTRFile,'The trial step is not good, rho = %1.4e\n',rho);
        fprintf(option.SubTRFile,'ared = %1.3e   pred = %1.3e\n',ared, pred);
        
% $$$         if (abs(pred) > 1.e5), keyboard, end
        % if the actual reduction is very small or there is no
        % reduction at all, shrink the TR size
        % and find the new acceptable point on the piesewise linear
        % path that already computed (BackTracking)
        
        trcountMax = ceil(log2(Delta/option.TolX));;%ceil(log2(Delta/norm(dirs(:,1))));
        if trcountMax == 0, trcountMax = 1; end
        trcount=1; lambda = 1;
        while (ared>=0) && (trcount <= trcountMax)
                
            %% Update the TR size Delta
             if (rho <= 0), Delta = gam0*min(Delta,norm(s));
             else Delta = min(gam1*Delta,gam1*norm(s)); end
            %Delta = gam1*norm(p);
            
            [s, dirs]=tradj(Delta, dirs); 
            
            pred = gc'*s + 0.5*s'*(feval(HessVect,xc,s));
             
            if (norm(s) < option.TolX), exitflag =2; histout(end,3) = norm(s);
                histout(end,4:5) = 0; OptimizationDisplay(histout,option); ...
                      
                    return; end
            
                if (pred > 0.0), pred, error('pred cannot be positive!'); end 
                 
                %fn = feval(f,xc+s);
            while 1,
                fn = feval(f,xc+s); 
                if ~isreal(fn)
                    s = s*0.8;
                else
                    break;
                end
            end
            if (abs(fn-histout(end,1)) <= option.TolF*(1+abs(histout(end,1)))), 
                exitflag = 1; xc = xc + s; histout(end,3) = norm(s);
                OptimizationDisplay(histout,option); return; end
        
            ared = fn-fc;
            
            %      Only compute a new pred if ared < 0 and there's some hope
            %      that rat > mu0
            if ared < 0
                pred=gc'*s + 0.5*s'*(feval(HessVect,xc,s));% + C.*s);
            end
            
            rho = ared/pred; trcount=trcount+1;
        end
        
        if ((rho > mu) & (ared < 0.0)) | (ared < 0.0)
            xc = xc + s; if (nargout > 4), XC(:,end+1) = xc; end
            if (rho<= beta) & (norm(s) < gam1*Delta),
                Deltao=Delta; Delta = gam1*Delta;
            elseif (rho>nu)
                if (rho > Laml)&&(norm(s) >= 0.8*Delta),
                    Deltao=Delta; Delta = min(gam2*Delta,Lamu);
                else
                    Deltao=Delta; Delta = min(max(Delta,gam2*norm(s)),Lamu);
                end
% $$$             elseif (rho>nu) & (norm(ps) >= 0.8*Delta), 
% $$$                 Deltao = Delta; Delta = min(gam2*Delta,Lamu);
            else, Deltao = Delta; end
        else
            rho, mu, norm(s), ared, Delta
            error(['Impossible, rho must be greater than mu after the ' ...
                   'while loop']); 
        end
                
    end
    histout(end,3:4) = [norm(s), Deltao]; histout(end,5) = CGIter;
    
    OptimizationDisplay(histout,option); % Display the Optimzation run

    if (histout(end,3) < option.TolX), exitflag = 2; return; end

    if (Delta < option.TolX), exitflag = 3; return; end

    it = it + 1; if (it > option.MaxIter), exitflag = 4; return; end

end
%
function OptimizationDisplay(histout,option)
it = size(histout,1);
fprintf(option.SubTRFile,'%d\t\t%2.3e\t%2.3e\t%2.3e\t%2.3e\t%d \n',...
    it-2,histout(end,1),histout(end,2),histout(end,3),histout(end,4),...
    histout(end,5));

%
% find the point of intersetion of the TR boundary and the piecewise
% linear path found in trcg
%
function [st, dirs] = tradj(trrad, dirs)
st=dirs(:,1); inside=1; itsl = size(dirs,2); kout = 0;
if norm(st) > trrad | itsl == 1
    st=st*trrad/norm(st);
    kout=1;
else
    for k=2:itsl
        if norm(st+dirs(:,k)) > trrad  & inside == 1
            kout=k;
            p=dirs(:,k); ac=p'*p; bc=2*(st'*p); cc=st'*st - trrad*trrad;
            alpha=(-bc + sqrt(bc*bc - 4*ac*cc))/(2*ac); st=st+alpha*p;
            inside=0; break;
        else
            st=st+dirs(:,k);
        end
    end
end
if (kout == 0), norm(st), trrad, 
  error('the TR radius cannot be bigger than norm(st)'); end
dirs = dirs(:,1:kout);
%-----------------------------------------------------------
%% Inexact Trust region Newton-CG 
function [x,directions,CGIter]  = NewtonCG(xc, f, gc, delta,pcv,HessVect)
%
% Solve the trust region problem with preconditioned conjugate-gradient
%
% C. T. Kelley, January 13, 1997
%
% Modified by Tan Bui, 2007, MIT
%
% This code comes with no guarantee or warranty of any kind.
% function [x, directions]
%                    = trcg(xc, f, gc, delta)
%
%
%
% Input:        xc=current point
%               b=right hand side
%           f = function, the calling sequence is
%                               [fun,grad]=f(x)
%           gc = current gradient
%                gc has usually been computed
%                before the call to dirdero
%           delta = TR radius
%           params = two dimensional vector to control iteration
%                params(1) = relative residual reduction factor
%                params(2) = max number of iterations
%           pcv, a routine to apply the preconditioner
%                if omitted, the identity is used.
%                The format for pcv is
%                       function px = pcv(x).
%
% Output:   x = trial step
%           directions = array of search directions TR radius reduction
%

%
% initialization
%
%%%----------------------------------------------------

n=length(xc); errtol = 1.e-3; maxiters = max(20,floor(n/2));
directions=zeros(n,1); NegCurve = 0; InsideTR = 1;
% the initial guess for Newton method is always zeros which means
% the initial guess is at the current point xc
x=zeros(n,1); b=-gc; r=b; 

% Check whether we should have preconditioner or not
if strcmpi(pcv,'')
    z=r;
else
    z = feval(pcv, xc,r);
end

rho=z'*r; tst=norm(r);
% terminate=errtol*norm(b);
terminate=min(1.e-3,norm(b))*norm(b); CGIter=1;

% Inexact Newton-CG loop with preconditioner
while((tst > terminate) & (CGIter <= maxiters) & (InsideTR == 1))
    
    if(CGIter==1), p = z; else
        beta=rho/rhoold; p = z + beta*p; end
    % approximate the product Hc*p
    %if (norm(p.*D)>50), keyboard, end
    w = feval(HessVect,xc,p);
    %w=H*w;
    alpha = p'*w;
    %
    % If alpha <=0 (search direction is in negative curvature direciton),
    % head to the TR boundary and return (return happens because
    % norm(x) will definitely be greater than hatdel)
    
    if(alpha <= 0)
        % Find alpha such that |x+alpha*p| = delta
        ac=p'*p; bc=2*(x'*p); cc=x'*x - delta*delta;
        alpha=(-bc + sqrt(bc*bc - 4*ac*cc))/(2*ac); NegCurve = 1;
        if (alpha < 0.0), error('Alpha must be positive!'); end
        %disp(' negative curvature')
        %delta
        InsideTR = 0;
    else
        % if alpha > 0 find the usual optimal step in the conjugate
        % direction
        alpha=rho/alpha;
        % If the new point is outside the TR, shrink it back
        if norm(x+alpha*p) > delta,
            ac=p'*p; bc=2*(x'*p); cc=x'*x - delta*delta;
            alpha=(-bc + sqrt(bc*bc - 4*ac*cc))/(2*ac);
            if (alpha < 0.0), error('Alpha must be positive!'); end
            InsideTR = 0;
        end

    end
    x=x+alpha*p; directions(:,CGIter)=alpha*p;
    r = r - alpha*w;
    % in CG, the residual is the negative of the gradient of the
    % approximate quadratic function
    tst=norm(r); rhoold=rho;
    if strcmpi(pcv,''), z=r; else z = feval(pcv,xc, r); end
    rho=z'*r; CGIter=CGIter+1;
end
% check the search direction
if (norm(x) > (delta+1.e2*eps))
    norm(x), delta
    error('norm(x) cannot be bigger than delta')
end
CGIter=CGIter-1;
