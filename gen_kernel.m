function[H,C,gamma] = gen_kernel(X,z)
% gen_kernel
% 
% first run chkbd.m or other program to generate
% m data points in X(:,1:m) with labels z(1:m)
% Then generate Gaussian kernel with parameters
% gamma and C as specified next.

gamma=0.03;  % constant in Gaussian kernel   exp(-gamma*||x-y||^2)

C=1e12;   % constant for the ``soft margin''

[~,m] = size(X);

K = eye(m);
for i = 1:m
    for j = 1:m
        K(i,j) = exp(-gamma*norm(X(:,i)-X(:,j))^2);
    end
end
H = diag(z)*K*diag(z); % scaled kernel

% Now, solve min 0.5*x'*H*x -e'*x | x'*z = 0; C*e >= x >= 0   (***)
%
% Here, z is a {-1,1} - vector !!!
% The variable x in (***) is not to be confused with the data matrix X

end
