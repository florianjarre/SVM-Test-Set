function beta = findbeta(x,H,z,C)
% by complementarity the parameter mu in the optimality conditions 
% coincides with the parameter beta used in the classification of the SVM

[n,~] = size(H);
c=ones(n,1);
g=H*x-c;
tolAS=1.0e-8; % use only components that are clearly inactive
INACT = x >= tolAS & x <= (1-tolAS)*C;
tmp = g.*z;
if sum(INACT) > 0
    beta = sum(tmp(INACT))/sum(INACT);
else
    warning('No reasonable value of beta found');
    sigma = zeros(size(g));
    sigma(x<tolAS) = 1;
    sigma(x>(1-tolAS)*C) = -1;
    sg = sigma.*g;
    sz = sigma.*z;
    beta = min(sg(sz==1));
end