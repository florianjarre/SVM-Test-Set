% testhalfmoon.m
% after solving the dual kernel problem determine the associated value
% beta and then TEST the regions divided by the svm using ntest points
% Input: H, c, z and the solution x from
%    min 0.5*x'*H*x -c'*x | z'*x = 0; C*e >= x >= 0.
% ( C is not used )

ntest = 10000; % number of test points 

[Xtest,ztest] = halfmoon(d,delta,rho,ntest);
% generate test data from the same distribution

beta = findbeta(x,H,z,C); 

[~,n] = size(H);
zR = ztest;
for ii = 1:ntest
    xtilde = Xtest(:,ii);
    kappatx = 0;
    for i = 1:n
      kappatx = kappatx + c(i)*x(i)*z(i)*exp(-gamma*norm(X(:,i)-xtilde)^2);
    end
    ztilde = sign(kappatx-beta);
    zR(ii) =  ztilde;
end

nh = ceil(ntest/2);
err1 = sum(ztest(1:nh)-zR(1:nh))/ntest;
err2 = sum(ztest(nh+1:ntest)-zR(nh+1:ntest))/ntest;
disp([err1,err2])

