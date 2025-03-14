function [x,outfile] = rsmo(H,z,C)
% Randomized SMO method for the problem
%    min 0.5*x'*H*x -c'*x | z'*x = 0; C*e >= x >= 0.
% with c=e in R^n. 
%
% Input:
% H n-by-n, positive definite, z in {-1,1}^n, C>0 in R
% 
% Output:
% final iterate x 

tolSMO = 1.0e-8;            % stopping tolerance 
dispkktnorm = 1;            % 1: print accuracy

% Initialization
[n,~] = size(H);
c = ones(n,1);
maxit = 10000*n;            % maximum number of iterations 
dH  = diag(H);


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


iter = 0;
happy = false;
violkkt = Inf;

while ~happy % begin OUTER ITERATIONS
    
    iter   = iter + 1;
    % determine i,j
    i = randi(n,1);
    j = randi(n,1);
    if i == j
        j = mod(i+1,n);
        if j == 0
            j = 1;
        end
    end % search direction e_i - z_i z_j e_j
    
    
    gij = g(i)-z(i)*z(j)*g(j);
    Hij = dH(i)+dH(j)-2*z(i)*z(j)*H(j,i); 
    lambdahat = -gij/Hij;
    if lambdahat > 0
        lambdahat = min(lambdahat,C-x(i));
    else
        lambdahat = max(lambdahat, -x(i));
    end
    if -z(i)*z(j)*lambdahat > 0
        if lambdahat > 0
            lambdahat = min(lambdahat,C-x(j));
        else
            lambdahat = max(lambdahat,x(j)-C);
        end
    else
        if lambdahat > 0
            lambdahat = min(lambdahat,  x(j));
        else
            lambdahat = max(lambdahat, -x(j));
        end
    end
    lambdajbar = lambdahat;
    
    % update x, g
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
    

    if mod(iter,10*n) == 0
        kkt_norm % compute  violkkt
    end
        
    
    
    % test for convergence
    if iter >= maxit || violkkt < tolSMO
        happy = true;
    end
    
end % end OUTER ITERATIONS

kkt_norm
if dispkktnorm
    disp('KKT violations: bounds, constraint, Lagrange gradient')
    disp([kviol1,kviol2,kviol3,kviol4])
end

outfile.iter = iter;
outfile.kviol1 = kviol1;
outfile.kviol2 = kviol2;
outfile.kviol3 = kviol3;
outfile.kviol4 = kviol4;
outfile.qnew = x.'*(H*(0.5*x)-c);
outfile.max_x = max(x);

end
