function[X,z] = halfmoon(d,delta,rho,n)
% halfmoon.m
% generate n sample points x^i \in R^d and labels z_i \in {\pm 1} 
% indicating whether x^i lies in some halfmoon shaped set or not.
%
% Input for example: d=3;delta=.25;rho=.25;n=100;
%                    rho > delta makes the separation easier
%
% First, define a half-moon-shape 
% S1 = { x \in R^d |  ||x||_2 <= 1, ||x+delta*e_1||_2 >= 1 }
% where delta \in (0,2) is a parameter.
% The support vector machine then is to decide whether a given "new"
% data point lies in S1 or not.
%
% Define S^- = S1 - rho*e_1 and S^+ = S1 + rho*e_1 and set 
% S2 = S^- union with S^+.
%
% The sample points x=x^i are then generated with some random distribution
% within the union of S1 and S2 and labels 
% z_i =  1 if x \in S1 and z_i = -1 if x \in S2.
% A point x in S1 is generated as follows:
% first x_1 is generated uniformly from [-delta/2,1],
% then a random normal vector y \in R^{d-1} is drawn and yn := y / ||y||_2,
% then x := [x_1; lambda*y] where lambda is drawn uniformly from
% [delta1, delta2] with 
% delta1 = sqrt(max{0, 1-(x_1+delta)^2}), delta2 = sqrt(1-x_1^2).
% Half of the points are then shifted to S2.

% first generate all n points in S1
x1 = (1+delta/2)*rand(1,n)-delta/2; % row vector in (-delta2,1)^n
y = randn(d-1,n); 
if d > 2
    yn = (sum(y.*y)).^(-0.5); % row vector with inverse column norms
else
    yn = 1./abs(y);
end
yn = y .* (ones(d-1,1)*yn); 
delta1 = max(0, 1-(x1+delta).^2).^.5;
delta2 = (1-x1.^2).^.5;
lambda = delta1 + rand(1,n).*(delta2-delta1);
X = [x1; yn .* (ones(d-1,1)*lambda)];

if max(sum(X.*X)) > 1
    error('generation in halfmoon failed 1')
end
tmpx = X; tmpx(1,:) = tmpx(1,:)+delta;
if min(sum(tmpx.*tmpx)) < 1
    error('generation in halfmoon failed 2')
end

% Now, move the last n/2 points to S2 ...
nh = ceil(n/2);
nq = ceil(n/4);
X(1,nh+1:nh+nq) = X(1,nh+1:nh+nq)-rho;
X(1,nh+nq+1:n)  = X(1,nh+nq+1:n) +rho;

% ... and assign the labels accordingly
z = ones(n,1);
z(nh+1:n) = -1;

end
