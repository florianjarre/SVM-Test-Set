function [x,outfile] = sweep(H,z,C)
% Active Set Sweeping method for the problem
%    min 0.5*x'*H*x -c'*x | z'*x = 0; C*e >= x >= 0
% with c=e in R^n. 
%
% Input:
% H n-by-n, positive definite, z in {-1,1}^n, C>0 in R
%
% Output:
% Some approximate solution x 
% and an estimate of its violation of the KKT conditions in kviol1 - kviol4
%
% experimental code 
% absolutely NO GUARANTEE for correctness, efficiency, or anything else
%
% Florian Jarre, Feb. 28, 2025

dispkktnorm = 1;           % 1: measure and print accuracy
refinelast = 1;            % 1: at the very end do one step of refinement
itrefine = 1;              % 1: at each iteration do a refinement step
tolAS  = 1.0e-14;          % tolerance for active set
regterm  = 1.0e-12*sum(abs(diag(H))); 
% regularization term to keep Cholesky & its updates well defined
% compensate with one or two steps of iterative refinement
          
% Initialization
[n,~] = size(H);
c = ones(n,1);
overallit = 0;
maxouterit = n;
wmsg1 = false; % no warning messge in sweep yet
wmsg2 = false; % no warning messge in up-cycle yet

% starting point x
g = -c;                    %                    gradient at x=0
s = -g + ((z.'*g)/n)*z;    % negative projected gradient at x=0
lambda = -g.'*s/(s.'*H*s); % minimizing q(lambda*s) 
lambda = min(lambda,C);    % this should rarely change lambda
x = lambda*s;              % intermediate point to define starting point
if max(x) > C
    lambda = min(C./s);     % use the fact that all s_i are positive
    x  = lambda*s;          % redefine the intermediate point
end
xtmp=x;

ACT1 = x < tolAS;          % determine the active set of x
ACT2 = x > (1-tolAS)*C;
ACT = ACT1 | ACT2;
INACT = ~ACT;
if sum(INACT) > 2000       % redefine sparser starting point 
    I1 = z>0;
    I2 = z<0;
    ncomp = min([sum(I1),sum(I2),1000]); %choose 2*ncomp entries from I1,I2
    g = H*x-c;
    s = -g + ((z.'*g)/n)*z;
    w = x+s;               % some mean value of projected gradient in [0,x]
    s = zeros(n,1);
    s(pilar(w,I1,ncomp)) = 1; % pick ncomp largest elemens of w(I1)
    s(pilar(w,I2,ncomp)) = 1;
    lambda = -g.'*s/(s.'*H*s);
    lambda = min(lambda,C); 
    x = lambda*s;
end

qnew = (0.5*x.'*H-c.')*x;  % computed function values may be low precision
sec_ref    = 0; % number      of second refinement steps (just to see)
no_sec_ref = 0; % number without second refinement steps (just to see)
outerit = 0;    % number of outer iterations
happy = false;  % overall termination criterion

nIN = sum(x >= tolAS & x <= C*(1-tolAS));
if nIN == 0 % first Newton step needs inactive components
    x=xtmp; % inefficient when C is small and n is large
end
clear xtmp

while ~happy % begin OUTER ITERATIOnS
outerit = outerit+1;
    
hitboundary = true;      % continue SWEEP as long as hitting the boundary 
sweepsteps = 0;          % number of steps in current sweep




while hitboundary        % Begin SWEEP 
sweepsteps = sweepsteps + 1;
if sweepsteps > n   
    error('sweep failed')
end

% one Newton step
if sweepsteps == 1   % refactor Cholesky only at the beginning of the sweep 
    ACT1 = x < tolAS;                      % determine the active set of x
    ACT2 = x > (1-tolAS)*C;
    ACT = ACT1 | ACT2;
    INACT = ~ACT;
    nIN = sum(INACT);
    AIN = H(INACT,INACT)+regterm*eye(nIN); % Solve A*s = b,  zI'*s = 0
    RA = chol(AIN);
end
b = c(INACT) - H(INACT,:)*x; 
zI = z(INACT);
xI = x(INACT);       % components strictly satisfying bound constraints
uv =  RA\(RA.'\ [b,zI]);                    % solve AIN*uv = [b,zI], 
if itrefine 
duv = RA\(RA.'\([b,zI]-H(INACT,INACT)*uv)); % iterative refinement
uv = uv+0.9*duv;     % reduced step length for the first correction step
if norm(duv,'fro') > 1e-3*norm(uv,'fro')
    duv = RA\(RA.'\([b,zI]-H(INACT,INACT)*uv)); % second refinement step
    uv = uv+duv;  
    sec_ref = sec_ref+1;
else
    no_sec_ref = no_sec_ref+1;
end
end
mu = (zI.'*uv(:,1))/(zI.'*uv(:,2));
s = uv(:,1)-mu*uv(:,2);
alpha = 1;


% line search for inactive constraints only
if min(xI+s) <= 0 || max(xI+s) >= C % full step is at or beyound boundary
    tmp1 = xI+s<0;
    tmp2 = xI+s>C;
    alpha1 = 1; 
    alpha2 = 1; 
    if sum(tmp1) > 0 
        [alpha1,~] = min(-xI(tmp1)./s(tmp1)); 
    end
    if sum(tmp2) > 0 
        [alpha2,~] = min((C-xI(tmp2))./s(tmp2));
    end
    alpha = min([1,alpha1,alpha2]);
    if alpha > 1+tolAS || alpha <= -tolAS
        error('programming error in Newton step length')
    end
else
    hitboundary = false;
end
sfull = zeros(n,1);
sfull(INACT) = s;      % needed for the update of the active set
xI = xI+alpha*s;       % alpha can be zero
xI = max(xI,0);        % in case of rounding errors
xI = min(xI,C);        % in case of rounding errors
x(INACT) = xI;
qold = qnew;
qnew = (0.5*x.'*H-c.')*x;
if qnew > qold + 1e-5*abs(qold) && ~wmsg1 
                              % allow for fairly large rounding error
    warning('increase of objective function in sweep')
    disp(' Possibly poorly conditioned problem' )
    wmsg1 = true;
end

% update active set
if hitboundary
    if alpha1 <= alpha2 
        tmp = abs(sfull)>0 & x <= tolAS*(abs(sfull)+1);
        [~,inew] = max(tmp);          % a new active index
        ACT1(inew) = true;
        tmp = abs(s)>0 & xI <= tolAS*(abs(s)+1); 
                                      % exactly the same indices
        [~,iiinew] = max(tmp);        % position ofthe new active index
    else
        tmp = abs(sfull)>0 & x >= (1-tolAS)*C;
        [~,inew] = max(tmp);          % a new active index
        ACT2(inew) = true;
        tmp = abs(s)>0 & xI >= (1-tolAS)*C; % exactly the same indices
        [~,iiinew] = max(tmp);        % position of the new active index
    end
    INACT(inew) = false;
    RA = choldelete(RA,iiinew);
    % the update of only one index of the active set possibly results in
    % a step of length zero in the next iteration if in fact two
    % constraints have turned active during the line search
end

if sum(INACT) == n
    hitboundary = false; % stop Newton iterations
end

end                      % End of SWEEP



overallit = overallit+sweepsteps;
upfail = false;   % UP-CYCLE until failing to add inactive variables
upit = 0;         % iterations in up-cycle



while ~upfail     % begin UP-CYCLE

% move away from active set  
% determine the active set of x (not of xJ)
upit = upit + 1;
ACT1 = x < tolAS;
ACT2 = x > (1-tolAS)*C;
ACT = ACT1 | ACT2;
INACT = ~ACT;
sigma = zeros(n,1);
sigma(ACT1) = 1;
sigma(ACT2) = -1; % s with sigma.*s >= 0 and z'*x=0 is a feasible direction 
g = H*x-c;
if sum(INACT) > 0
    muu = g(INACT).'*z(INACT)/sum(INACT);
else
    muu = (z.'*g)/n;
end
gtilde = g -muu*z;
stilde = -gtilde;
stilde(ACT1) = max(0,stilde(ACT1));
stilde(ACT2) = min(0,stilde(ACT2));
stilde(INACT) = 0;
II = z.*stilde > 0;
JJ = z.*stilde < 0;
s = zeros(size(stilde));
if sum(II) > 0 && sum(JJ) > 0  % Case 1.
    v1 = z(II).'* stilde(II);
    v2 = z(JJ).'* stilde(JJ);
    s(II) = -v2*stilde(II);
    s(JJ) =  v1*stilde(JJ);
elseif sum(II)+sum(JJ) > 0     % Case 2.
    [~,p1] = max(abs(stilde)); % p1 in II or JJ since stilde(INACT) = 0
    tmps = zeros(size(stilde)); 
    tmpz = (-sigma(p1)*z(p1))*z.*sigma >= 0; % indices that complement p1
    tmps(tmpz) = (-sigma(p1)*z(p1))*z(tmpz).*gtilde(tmpz);
    tmps(~tmpz)   = Inf;
    [~,p2] = min(tmps);
    s(p1) =   sign(stilde(p1));
    s(p2) =  -sign(stilde(p1))*z(p1)*z(p2);
else
    upfail = true;
end

if g.'*s >= 0
    upfail = true;
    lambda = 0;
else
    lambda = -g.'*s/(s.'*H*s);
    tmp = s>0;
    if sum(tmp) > 0
        lambda = min([lambda;(C-x(tmp))./s(tmp)]);
    end
    tmp = s<0;
    if sum(tmp) > 0
        lambda = min([lambda;-x(tmp)./s(tmp)]);
    end
end
x = x + lambda*s;    
x = max(x,0);
x = min(x,C);
ACT = x < tolAS | x > (1-tolAS)*C;
qold = qnew;
qnew = (0.5*x.'*H-c.')*x;
if qnew > qold + 1e-5*abs(qold) && ~wmsg2
                           % allow for fairly large rounding error
    warning('increase of objective function in up-cycle')
    disp(' Possibly poorly conditioned problem' )
    wmsg2 = true;
end

if sum(ACT) >= n
    upfail = true; % not a fail, but finish the up-cycle
    warning('all indices active in up-cycle') % message for debugging
end
if upit > n
    upfail = true; % possible cycling
    warning('more than n steps in up-cycle') % message for debugging
end

end % end UP-CYCLE

if upit == 1
    happy = true; % after a sweep the first step of the up-cycle failed
end
if outerit >= maxouterit
    happy = true;
end
overallit = overallit+upit;

end % end OUTER ITERATIONS

kkt_norm
if dispkktnorm
    disp('KKT violations: bounds, constraint, Lagrange gradient')
    disp([kviol1,kviol2,kviol3,kviol4])
end





% In case of rounding errors use cheap correction to make x (more) feasible
ACT1 = x < tolAS;        % determine the active set of x
ACT2 = x > (1-tolAS)*C;
x(ACT1) = 0;
x(ACT2) = C;
ACT = ACT1 | ACT2;
INACT = ~ACT;
nIN = sum(INACT);
if nIN > 0
    muu = x.'*z/nIN;
    x(INACT) = x(INACT) - muu*z(INACT);
    x = max(x,0); % in case some inactive components got active
    x = min(C,x); % -- this might result in violating again z'*x = 0
    INACT = x >= tolAS & x <= (1-tolAS)*C;
    nIN = sum(INACT);
end

% Do a last Newton step to reduce the final error (copied from above)
if refinelast %-----------irefine (mostly copied from above)
if nIN > 0            % TYPICAL CASE: Newton with inactive coordinates
AIN = H(INACT,INACT); % Solve A*s = b,  zI'*s = 0
tmpt = 1000*eps*sum(diag(AIN));
RA = chol(AIN+tmpt*eye(nIN)); % refactor for higher accuracy
b = c(INACT) - H(INACT,:)*x; 
zI = z(INACT);
xI = x(INACT);       % starting point strictly satisfying bound constraints
uv = RA\(RA.'\[b,zI]);       % solve AIN*uv = [b,zI], 
duv = RA\(RA.'\([b,zI]-H(INACT,INACT)*uv)); % iterative refinement
uv = uv+0.9*duv;     % reduced step length for the first correction step
if norm(duv,'fro') > 1e-3*norm(uv,'fro')
    duv = RA\(RA.'\([b,zI]-H(INACT,INACT)*uv)); % iterative refinement
    uv = uv+duv;
    sec_ref = sec_ref+1;
else
    no_sec_ref = no_sec_ref+1;
end

mu = (zI.'*uv(:,1))/(zI.'*uv(:,2));
s = uv(:,1)-mu*uv(:,2);
alpha = 1;

% line search for inactive constraints only
if min(xI+s) <= 0 || max(xI+s) >= C % full step is at or beyound boundary
    tmp1 = xI+s<0;
    tmp2 = xI+s>C;
    alpha1 = 1; 
    alpha2 = 1; 
    if sum(tmp1) > 0 
        [alpha1,~] = min(-xI(tmp1)./s(tmp1));
    end
    if sum(tmp2) > 0 
        [alpha2,~] = min((C-xI(tmp2))./s(tmp2));
    end
    alpha = min([1,alpha1,alpha2]);
    if alpha > 1+tolAS || alpha <= -tolAS
        error('programming error in final Newton step length')
    end
end
x(INACT) = xI+alpha*s;  % alpha can be zero
x = max(x,0);           % in case of rounding errors
x = min(x,C);           % in case of rounding errors
qnew = (0.5*x.'*H-c.')*x;
end % of TYPICAL CASE: Newton with inactive coordinates
end %-----------irefine


kkt_norm
if dispkktnorm
    disp('KKT violations: bounds, constraint, Lagrange gradient')
    disp([kviol1,kviol2,kviol3,kviol4])
end

outfile.outerit = outerit;      % number of cycles
outfile.overallit = overallit;  % number of up-steps + number of down-steps 
outfile.kviol1 = kviol1;        % violation of x >= 0
outfile.kviol2 = kviol2;        % violation of x <= C
outfile.kviol3 = kviol3;        % relative violation of z'*x = 0
outfile.kviol4 = kviol4;        % relative norm of projected gradient
outfile.no_sec_ref = no_sec_ref;% iterations without second refinement step
outfile.sec_ref = sec_ref;      % iterations with    second refinement step
outfile.qnew = qnew;            % final function value
outfile.max_x = max(x);         % largest component of x

end
