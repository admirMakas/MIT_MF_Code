function [xc,fc,exitflag,histout,XC] = STIRNCG(xc,f,low,up,option)
%  STIRNCG  Copyright and License information for STIRNCG  
%      STIRNCG is the subspace trust region interior reflective Newton-CG method
% 
%      The following information is a copy of the License file in the STIRNCG
%      distribution.
%   
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   
%    STIRNCG, Copyright (c) 2007 by Tan Bui
%    All Rights Reserved. Massachusetts Institute of Technology
%    
%    STIRNCG License:
%    
%        Your use or distribution of STIRNCG or any modified version of
%        STIRNCG implies that you agree to this License.
%    
%        THIS MATERIAL IS PROVIDED AS IS, WITH ABSOLUTELY NO WARRANTY
%        EXPRESSED OR IMPLIED.  ANY USE IS AT YOUR OWN RISK.
%    
%        Permission is hereby granted to use or copy this program, provided
%        that the Copyright, this License, and the Availability of the original
%        version is retained on all copies.  User documentation of any code that
%        uses STIRNCG or any modified version of STIRNCG code must cite the
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
%    A Subspace, Interior, and Conjugate Gradient Method for Large-Scale 
%    Bound-Constrained Minimization Problems. SIAM Journal on Scientific
%    Computing, Volume 21 ,  Issue 1  (Aug.-Sept. 1999), Pages: 1 - 23
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
global it OutputFile

% hardcoded all the parameters
mu = 1.e-1; beta = 0.25; nu = 0.75; gam0 = 0.0625; gam1 = 0.5; gam2 = 2.0;
Laml = 1; maxits = 5; BackTrack = 0.7; 
n = length(low);
Lamu = max(sqrt(sum(min((up-low).^2,1000))),1);

if option.Scaled == 1
    sig = 1;
elseif option.Scaled == 0
    sig = 0.99995;
else
    option.Scaled
    error('unsupported Scaled')
end


if (nargin<5),
    option.TolF  = 1.e-12; option.TolX  = 1.e-10;
    option.TolGrad = 1.e-8; option.MaxIter = 100;
    option.Scaled = 1; option.OptMethod = 'STIRNCG';
    option.AffineScaling = 'huuIa'; option.ReflectionFlag = 1; 
    option.pcv = @Precond;
end

try
  OutputFile = option.OutputFile;
catch
  OutputFile = 1;
end

ReflectionFlag = option.ReflectionFlag;

% put initial iterate in feasible set
if norm(xc - proj(xc,up,low,option)) > 0
    fprintf(OutputFile,' initial iterate not feasibile, project it back to the feasible region \n');
    xc=proj(xc,up,low,option);
end

if sum(find((low-up)>0.0))>0, sum(find((low-up)>0.0))
    error('The Lower bound cannot be bigger than the upper bound'); end

it = 1; histout = zeros(1,5); XC = []; 
if (nargout > 4), XC(:,end+1) = xc; end

fprintf(OutputFile,'\n');
fprintf(OutputFile,'Iter \t\t  F \t\tScaled_g\t Step\t\t  TR size\t PosDef\n');

% The main loop of the iterative optimization
while (1)

    % compute the initial function and gradient
    [fc,gc]=feval(f,xc);
    histout(end+1,1) = fc;
    % Compute the Affine scaling
    [J,C,V,D] = AffineScaling(xc,gc,low,up,option);
    histout(end,2) = norm(gc.*V);
 
    % initial TR radius 
    %if (it == 1), Delta = 0.5e0*Lamu; Deltao = Delta; end
    if (it == 1), Delta = min(Lamu,histout(end,2)/5); Lamu = Delta; Deltao = Delta; end
    
    if ((it == 1)&(abs(histout(end,1)) < option.TolF)), histout(end,5)=0;
        OptimizationDisplay(histout, it-1); exitflag = 1; return; end
    % if scaling gradient is small, exit
    if (histout(end,2) < option.TolGrad), histout(end,5) = histout(end-1,5);
        OptimizationDisplay(histout, it-1); exitflag = 5; return; end

    % Inexact Newton CG for the scaled TR subproblem to determine
    % the inexact Newton direction or the negative direction
    [p,d,posdef,kcg] = MPCG(xc,f,gc,D,J,option);
    histout(end,5) = posdef;
    
    % Determine the subspace
    [S] = Subspace(p,d,posdef,kcg,D,gc,xc, option.HessVect);
    
    % Determine the coeffs of the subspace problem
    [a,b,c] = SubCoeffs(S,gc,xc,D,J, option.HessVect);
    
    % Solve the subspace trust region method
    SubspaceTRSolve(a,b,c,posdef,Delta, option.SubTRFile);
    
    % Solution
    load SubTRSolution
    
    ps = S*z; p = ps.*D; %clear z S a b c alpha d posdef
    
    
    % find the reflected direction 
    [pr,rflag]  = reflection(xc,p,up,low,gc.*V, ReflectionFlag);
    
    % Project both p and pr back to the interior
     if option.Scaled == 1
         p = sig*(proj(xc+p,up,low,option)-xc);
     elseif option.Scaled == 0
         xp = proj(xc+p,up,low,option);
         sigk = max(sig,1-norm(xp-xc));
         p = sigk*(xp-xc); 
     end
    
    psip = gc.'*p + 0.5*p'*(feval(option.HessVect,xc,p));% + C.*p);
    
    if (rflag == 1)
        fprintf(OutputFile,'Reflection\n');
        pr = sig*(proj(xc+pr,up,low,option)-xc);
        psipr = gc'*pr + 0.5*pr'*(feval(option.HessVect,xc,pr));% + C.*pr);
       
        if (psip < psipr), s = p; pred = psip; 
        else, s = pr; pred = psipr; end
     
    else, s = p; pred = psip; end
    
    if (norm(s) < option.TolX), exitflag = 2; histout(end,3) = norm(s);
        histout(end,4:5) = [Delta,posdef]; OptimizationDisplay(histout,it-1); return; end
% $$$     try
% $$$      if (pred > 0.0), pred, error('pred cannot be positive!'); end
% $$$     catch, keyboard, end
    
    % actual reduction versus predicted reduction
    fn = feval(f,xc+s); 
    if (abs(fn-histout(end,1)) <= option.TolF*(1+abs(histout(end,1)))), 
            exitflag = 1; xc = xc + s; histout(end,3) = norm(s);
       histout(end,4:5) = [Delta,posdef]; OptimizationDisplay(histout,it-1); return; end
        
    ared = fn-fc;% + 0.5*s'*(C.*s);
    rho = ared/pred;
   
    % If the actual reduction is reasonable, take the new point
    % and adjust the TR radius and move on
    if ((rho > mu) & (ared < 0.0)) | (ared < 0.0)
        xc = xc + s; if (nargout > 4), XC(:,end+1) = xc; end
        if (rho<= beta) & (norm(ps) < gam1*Delta),
            Deltao=Delta; Delta = gam1*Delta;
        elseif rho > nu,
            if (rho > Laml)&&(norm(ps) >= 0.8*Delta),
                Deltao=Delta; Delta = min(gam2*Delta,Lamu);
            else
               Deltao=Delta; Delta = min(max(Delta,gam2*norm(ps)),Lamu);
            end
% $$$         elseif (rho>nu) & (norm(ps) >= 0.8*Delta), 
% $$$             Deltao = Delta; Delta = min(gam2*Delta,Lamu);
        else, Deltao = Delta; end

    else
        fprintf(OutputFile,'The trial step is not good, rho = %1.4e\n',full(rho));
        fprintf(OutputFile,'ared = %1.3e   pred = %1.3e\n',full(ared), full(pred));
        %if (abs(pred) > 1.e6), keyboard, end
        % if the actual reduction is very small or there is no
        % reduction at all, shrink the TR size
        % and find the new acceptable point on the piesewise linear
        % path that already computed (BackTracking)
        
        trcountMax = ceil(log2(Delta/option.TolX));;%ceil(log2(Delta/norm(dirs(:,1))));
        if trcountMax == 0, trcountMax = 1; end
        trcount=1; lambda = 1;
        while (ared>=0) && (trcount <= trcountMax)
                
            %% Update the TR size Delta
             if (rho <= 0), Delta = gam0*min(Delta,norm(ps));
             else Delta = min(gam1*Delta,gam1*norm(ps)); end
            %Delta = gam1*norm(ps);
            
            % Solve the subspace trust region method
            SubspaceTRSolve(a,b,c,posdef,Delta, option.SubTRFile);
            
            % Solution
            load SubTRSolution
            
            ps = S*z;  
            
            if option.Scaled == 1, p = ps.*D; else, p = ps; end
            
            [pr,rflag]  = reflection(xc,p,up,low,gc.*V, ReflectionFlag);
            
            % Project both p and pr back to the interior
            if option.Scaled == 1
                p = sig*(proj(xc+p,up,low,option)-xc);
            elseif option.Scaled == 0
                xp = proj(xc+p,up,low,option);
                sigk = max(sig,1-norm(xp-xc));
                p = sigk*(xp-xc); 
            end
   
            psip = gc'*p + 0.5*p'*(feval(option.HessVect,xc,p));% +
                                                          % C.*p);
             
             if (rflag == 1)
                 fprintf(OutputFile,'Reflection\n');
                 pr = sig*(proj(xc+pr,up,low,option)-xc);
                 psipr = gc'*pr + 0.5*pr'*(feval(option.HessVect,xc,pr));% + C.*pr);
                 
                 if (psip < psipr), s = p; pred = psip;
                 else, s = pr; pred = psipr; end
 
             else, s = p; pred = psip; end
             % s = p;
             if (norm(s) < option.TolX), exitflag =2; histout(end,3) = norm(s);
                  histout(end,4:5) = [Delta,posdef]; OptimizationDisplay(histout,it-1); ...
                      
                    return;
             end
            
% $$$                  try
% $$$                      if (pred > 0.0), pred, error('pred cannot be positive!'); end
% $$$                  catch, keyboard, end   
                 
            fn = feval(f,xc+s);
            if (abs(fn-histout(end,1)) <= option.TolF*(1+abs(histout(end,1)))), 
                exitflag = 1; xc = xc + s; histout(end,3) = norm(s);
                histout(end,4:5) = [Delta,posdef]; OptimizationDisplay(histout,it-1); return; end
        
            ared = fn-fc;% + 0.5*s'*(C.*s);
            
            %      Only compute a new pred if ared < 0 and there's some hope
            %      that rat > mu0
            if ared < 0
                pred=gc'*s + 0.5*s'*(feval(option.HessVect,xc,s));% + C.*s);
            end
            
            rho = ared/pred; trcount=trcount+1;
        end
        
        if ((rho > mu) & (ared < 0.0)) | (ared < 0.0)
            xc = proj(xc + s,up,low,option); if (nargout > 4), XC(:,end+1) = xc; end
            if (rho<= beta) & (norm(ps) < gam1*Delta),
                Deltao=Delta; Delta = gam1*Delta;
            elseif (rho>nu)
                if (rho > Laml)&&(norm(ps) >= 0.8*Delta),
                    Deltao=Delta; Delta = min(gam2*Delta,Lamu);
                else
                    Deltao=Delta; Delta = min(max(Delta,gam2*norm(ps)),Lamu);
                end
            
            else, Deltao = Delta; end
        else
          exitflag = 6;
          break;
          rho, mu, norm(s), ared, Delta
          error(['Impossible, rho must be greater than mu after the ' ...
                 'while loop']); 
            
        end
                    
    end

    histout(end,3:4) = [norm(ps), Deltao]; histout(end,5) = posdef;
    
    OptimizationDisplay(histout,it-1); % Display the Optimzation run

    if (histout(end,3) < option.TolX), exitflag = 2; return; end

    if (Delta < option.TolX), exitflag = 3; return; end

    it = it + 1; if (it > option.MaxIter), exitflag = 4; return; end

end
%
function OptimizationDisplay(histout, it)
global OutputFile
fprintf(OutputFile,'%d\t\t%2.3e\t%2.3e\t%2.3e\t%2.3e\t%d\n',...
    it-1,histout(end,1),histout(end,2),histout(end,3),histout(end,4),...
    histout(end,5));
%------------------------------------------------------------
function [J,C,V,D] = AffineScaling(xc,gc,low,up,option)
n = length(gc); Index = 1:n';
V = zeros(n,1); I = find(gc<0.0); V(I) = abs(up(I) - xc(I));
I = setdiff(Index,I); V(I) = abs(xc(I) - low(I));
if (~isempty(find(V<(-eps)))), V, error('V < 0.0, check!'), end
Iinf = isinf(V);
V(Iinf) = 1;
C = 0;

% Compute the J and C matrix
if strcmpi(option.AffineScaling,'CL')
    % Coleman-Li scaling
    D = sqrt(V); J = sign(gc).*gc;  
    if option.Scaled == 0,
          C = J./V; 
    end
elseif strcmpi(option.AffineScaling,'CLModified')
    % Modified CL scaling
    D = V; J = 2*sign(gc).*D.*gc; V = V.^2; %C = J./V; 
elseif strcmpi(option.AffineScaling,'Tan')
    % My scaling
    D = V; J = sign(gc).*gc.*V; %C = J./V.^2;
elseif strcmpi(option.AffineScaling,'HUUIa')
    % HUU scaling Ia
    p = 2; D = ones(size(gc)); J = zeros(size(gc));
    I = union(find(abs(gc)<min(xc-low,up-xc).^p),...
        find(min(xc-low,up-xc)<abs(gc).^p));
    Dcl = sqrt(V); D(I) = Dcl(I); J(I) = abs(gc(I));
    V = D.^2; %C = J./V;
elseif strcmpi(option.AffineScaling,'HUUIb')
    % HUU scaling Ib
    p = 2; D = ones(size(gc)); J = zeros(size(gc));
    I = union(find(abs(gc)<min(xc-low,up-xc).^p),...
        find(min(xc-low,up-xc)<abs(gc).^p));
    Dcl = V; D(I) = Dcl(I); J(I) = abs(gc(I)).*D(I);
    V = D; %C = J./V.^2;   
elseif strcmpi(option.AffineScaling,'HUUIIa')
    % HUU scaling Ia
    p = 2; J = zeros(size(gc));
    I = union(find(abs(gc)<min(xc-low,up-xc).^p),...
        find(min(xc-low,up-xc)<abs(gc).^p));
    D = sqrt(V); J(I) = abs(gc(I)); %C = J./V;
elseif strcmpi(option.AffineScaling,'HUUIIb')
    % HUU scaling Ib
    p = 2; J = zeros(size(gc)); C = J;
    I = union(find(abs(gc)<min(xc-low,up-xc).^p),...
        find(min(xc-low,up-xc)<abs(gc).^p));
    D = V; C(I) =  abs(gc(I)); J(I) = C(I).*D(I); 
    if option.Scaled == 0,
         C = C./D;
    end
elseif strcmpi(option.AffineScaling,'HUU')
    % HUU scaling 
    p = 2; J = zeros(size(gc)); E = J;
    I = union(find(abs(gc)<min(xc-low,up-xc).^p),...
        find(min(xc-low,up-xc)<abs(gc).^p));
    D = V; E(I) = abs(gc(I)); W = 1./(D+E);
    D = W.*D;
    J(I) = abs(gc(I)).*D(I).*W(I); %C = J./V.^2;  
elseif strcmpi(option.AffineScaling,'Original')
    % No scaling
    D = ones(size(V)); V = D; J = zeros(size(V)); %C = J;
else, error('unsupported affine scaling'), end
D(Iinf) = 1; V(Iinf) = 1; J(Iinf) = 0;
%-------------------------------------------------------------    
% projection and scaling onto the strictly feasible set
function px = projs(xc,p,kku,kkl)
% $$$ epsilon = 1.e-12;
% $$$ px=min(kku-epsilon,x);
% $$$ px=max(kkl+epsilon,px);
x = xc+p;
I = find((x-kku) > (2*eps)); a = min(abs(kku(I)-xc(I))./p(I));
I = find((kkl-x) > (2*eps)); b = max(abs(xc(I)-kkl(I))./p(I));

if ~isempty(a)
    if (a<0.0), a, error('a must be positive'); end
    if ~isempty(b)
        if (b>0.0), b, error('b must be negative'); end
        c = min(a,-b);
    else
        c = min(a);
    end
else
    if ~isempty(b)
        if (b>0.0), b, error('b must be negative'); end
        c = min(-b);
    else
        c = 1;
    end
end

if (c>1), c, error('c must be less than 1'); end

if (c > 1.e-8)
    px = p*c;
else
    px = proj(x,kku,kkl)-xc;
end

%------------------------------------------------------------------
function [p,rflag] = reflection(xc,p,up,low,gc, ReflectionFlag)

xc = xc + p; rflag = 0;
I1 = find((up-xc) < 0.0); I2 = find((low-xc) > 0.0);
emptyflag = 0;

% single reflection
if (ReflectionFlag == 1)
    if (~isempty(I1))
        if (~isempty(I2))
            %         y = [(up(I1(1))-xc(I1(1))),(low(I2(1))-xc(I2(1)))];
            %[I,It] = min([I1(1),I2(1)]); %y = y(It);
            I = sort([I1;I2],'ascend');
        else
            I = I1; %y = (up(I1(1))-xc(I1(1)));
        end
    else
        if (~isempty(I2))
            I = I2;%min(I2); %y = (low(I2(1))-xc(I2(1)));
        else
            emptyflag = 1;
        end
    end
    if (emptyflag == 0),
        nI = length(I); 
        for ii = 1:nI
            if (gc(I(ii))*p(I(ii)) > 0.0)
                p(I(ii)) = -p(I(ii)); rflag = 1;
                break; 
            end
        end
    end
elseif (ReflectionFlag == 2)
    % Reflection all the descent direction
    I = unique([I1;I2]);
    if (~isempty(I))
        n = length(I);
        for in =1:n
            if (gc(I(in))*p(I(in)) > 0)
                p(I(in)) = -p(I(in)); rflag = 1;
            end
        end
    end
    
elseif (ReflectionFlag == 0)
    % Do nothing 
else
    ReflectionFlag
    error('unsupported ReflectionFlag');
end

%----------------------------------------------------------------------
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
% The modified  preconditioned conjugate gradient method for scaled
% problem 
function [p,d,posdef,k] = MPCG(xc,f,gc,D,J, option)

pcv = option.pcv;
HessVect = option.HessVect;
    
n = length(gc);  kmax = max(ceil(n/2),2); posdef = 1;
k = 0; p = zeros(n,1); r = -gc.*D; epsilon = sqrt(eps);
[P,Ri] = feval(pcv,xc); normr = 1;
eta = min(1.e-1,norm(r))*norm(Ri.*r);
rho = 0;
while k<kmax
   z = r./P; rhoold = rho; rho = r.'*z;
   k = k + 1;
   if k == 1,
       d = z;
   else
       beta = rho/rhoold;
       d = z + beta*d;
   end
   w = feval(HessVect,xc,d.*D).*D + J.*d;
   gam = d.'*w; temp = d.'*(P.*d);
   %if (gam < epsilon*temp), posdef = 0; return;
   if (gam <= 0.0), posdef = 0; return;
   else
       alpha = rho/gam;
       p = p + alpha*d;
       r = r-alpha*w;
   end
   if ((normr*norm(r)) <= eta)
       d = zeros(n,1); return; 
   end
       
end
d = zeros(n,1);

%----------------------------------------------------------
% Determine the subspace
function [S] = Subspace(p,d,posdef,k,D,gc,xc, HessVect)

tau = sqrt(eps);

if posdef == 1,
    if k == 1,
        S =  D.*sign(gc);%p;%
    else
      S = [D.*gc,p];
        %[S,R] = qr([D.*gc,p],0);
% $$$         S = D.*gc/norm(D.*gc);
% $$$         [p,nnv] = MGS(S,p);
% $$$         if (nnv ~= 0)
% $$$           S = [S,p];
% $$$         end
    end
elseif posdef == 0,
    w = d.*D;
    temp = D.*(D.*sign(gc));
    temp1 = temp.'*(feval(HessVect,xc,temp));
    temp2 = (norm(D.*(D.*gc))/norm(w))^2*(w.'*feval(HessVect,xc,w));
    if temp1 < (tau*temp2),
        S = D.*sign(gc);
    else
        if k == 1
            S = D.*sign(gc);%d;
        else
          S = [D.*sign(gc),d];
            %[S,R] = qr([D.*sign(gc),d],0);
% $$$             S = D.*sign(gc)/norm(D.*sign(gc));
% $$$             [d,nnv] = MGS(S,d);
% $$$             if (nnv ~= 0)
% $$$               S = [S,d];
% $$$             end
        end
    end
else
    posdef
    error('posdef must be 1 or 0')
end

if (size(S,2) == 2)
  [q,r]=qr(S);
  if (abs(r(2,2)) < tau*abs(r(1,2)))
    S = S(:,1);
  end
end

%------------------------------------------------------------------
% Determine the coeffs of the subspace problem
function [a,b,c] = SubCoeffs(S,gc,xc,D,J, HessVect)

n = size(S,2); gh = D.*gc;

if n == 1,
    temp = S.*D;
    a = temp.'*(feval(HessVect,xc,temp))+S.'*(J.*S);
    b = gh.'*S;  c = S.'*S;
elseif n == 2,
   a = zeros(3,1); b = zeros(2,1); c = zeros(3,1);
   
   temp1 = S(:,1).*D; temp1n = feval(HessVect,xc,temp1);
   a(1) = temp1.'*(temp1n)+S(:,1).'*(J.*S(:,1));
   temp2 = S(:,2).*D;
   a(3) = temp2.'*(feval(HessVect,xc,temp2))+S(:,2).'*(J.*S(:,2));
   a(2) = temp2.'*(temp1n)+S(:,2).'*(J.*S(:,1));
   
   b = S.'*gh; 
   
   c(1) = S(:,1).'*S(:,1); c(2) = S(:,1).'*S(:,2); 
   c(3) = S(:,2).'*S(:,2);
else
    n
    error('n must be 1 or 2')
end
