function [x,outfile] = gsmo(H,z,C)
% Greedy SMO method for the problem
%    min 0.5*x'*H*x -c'*x | z'*x = 0; C*e >= x >= 0.
% with c=e in R^n. 
%
% Input:
% H n-by-n, positive definite, z in {-1,1}^n, C>0 in R
% 
% Output:
% final iterate x 

tolSMO = 1.0e-10;           % stopping tolerance 
tolAS  = 1.0e-14;           % tolerance for active set
dispkktnorm = 1;            % 1: print accuracy

% Initialization
[n,~] = size(H);
c = ones(n,1);
dH  = diag(H);
maxit = 1000*n;             % maximum number of iterations 


% starting point x
g  = -c;                    %                    gradient at x=0
s  = -g + ((z.'*g)/n)*z;    % negative projected gradient at x=0
lambda = -g.'*s/(s.'*H*s);  % minimizing q(lambda*s) 
lambda = min(lambda,C);     % this should rarely change lambda
x  = lambda*s;              % starting point
if max(x) > C
    lambda = min(C./s);     % use the fact that all s_i are positive
    x  = lambda*s;          % redefine the starting point
end
g  = H*x-c;                 % gradient at x
ACT1 = x <= tolAS;
ACT2 = x >= (1-tolAS)*C;
INACT = ~(ACT1|ACT2);
nIN = sum(INACT);
if nIN > 0
    mu = g(INACT).'*z(INACT)/nIN;
    gt = g-mu*z; % KKT implies gt(ACT1) >= 0,  gt(ACT2) <= 0
else
    gt = g-(g.'*z/n)*z;
end




iter = 0;
lambdabar = zeros(n,1);
happy = false;
violkkt = Inf;

while ~happy % begin OUTER ITERATIONS
    
    iter   = iter + 1;
    % determine i
    st = -gt;      % stilde, descent direction not violating bounds
    if sum(ACT1) > 0
        st(ACT1) = max(0,st(ACT1));
    end
    if sum(ACT2) > 0
        st(ACT2) = min(0,st(ACT2));
    end
    [~,i] = max(abs(st));
    if st(i) < 0  % st(i) = 0 only when x is optimal
        lambdaibar = -x(i); % (negative)
    else
        lambdaibar = C-x(i); 
    end % feasible moves: x --> x+lambda e_i with lambda in <0,lambdaibar>
    
    % determine j
    if z(i) > 0
        gij = gt(i)-z.*gt;
    else
        gij = gt(i)+z.*gt;
    end
    tmp = gt(i).*gij;
    tmp(i) = -1;             % do not use index i
    candj = tmp > 0;         % candidates for choosing j
    Hij = dH(i)+dH-(2*z(i))*z.*H(:,i); % todo: eliminate n multiplications
    Hij(i) = dH(i);          % just to avoid NaN when forming lambdahat
    lambdahat = -gij./Hij;
    ind1 = -z(i)*z > 0;      % x(j) --> x(j) - z(i)*z(j)*lambdabar(j)
    ind2 = -z(i)*z < 0;
    tmp1 = min(lambdahat(ind1), C-x(ind1));
    tmp1 = max(tmp1,             -x(ind1));    
    lambdabar(ind1) = tmp1;
    tmp2 = min(lambdahat(ind2),   x(ind2));
    tmp2 = max(tmp2,            x(ind2)-C);    
    lambdabar(ind2) = tmp2;
    lambdabar = sign(lambdabar).*min(abs(lambdaibar),abs(lambdabar));
    delta = lambdabar.*(gij+(0.5*lambdabar).*Hij);
    delta(~candj) = 1; % only negative entries of delta lead to descent
    [deltamin,j] = min(delta);
    lambdajbar = lambdabar(j);
    if lambdaibar*lambdajbar < 0
        error('errror in step length comutation') % for debugging
    end
    
    % update x, g, gt
    x(i) = x(i) + lambdajbar;
    if min(x) < 0 || max(x) > C
        error('xi out of bound') % for debugging
    end
    x(j) = x(j) - z(i)*z(j)*lambdajbar;
    if min(x) < 0 || max(x) > C
        error('xj out of bound') % for debugging
    end
    if z(i)*z(j) > 0
        g = g + lambdajbar*(H(:,i)-H(:,j)); 
    else
        g = g + lambdajbar*(H(:,i)+H(:,j));
    end
    ACT1 = x <= tolAS;
    ACT2 = x >= (1-tolAS)*C;
    INACT = ~(ACT1|ACT2);
    nIN = sum(INACT);
    if nIN > 0
        mu = g(INACT).'*z(INACT)/nIN;
        gt = g-mu*z; % KKT implies gt(ACT1) >= 0,  gt(ACT2) <= 0
    else
        gt = g-(g.'*z/n)*z;
    end

    if mod(iter,10*n) == 0
        kkt_norm % compute  violkkt
    end
        
    
    
    % test for convergence
    if iter >= maxit || violkkt < tolSMO
        happy = true;
    end
    if deltamin > 0
        error('non-descent-step computed') % for debugging
    end
    
end % end OUTER ITERATIONS

kkt_norm
if dispkktnorm
    disp([kviol1,kviol2,kviol3,kviol4])
end

outfile.iter = iter;              % number of iterations
outfile.kviol1 = kviol1;          % violation of x >= 0
outfile.kviol2 = kviol2;          % violation of x <= C
outfile.kviol3 = kviol3;          % relative violation of z'*x = 0
outfile.kviol4 = kviol4;          % relative norm of projected gradient
outfile.qnew = x.'*(H*(0.5*x)-c); % final function value
outfile.max_x = max(x);

end
