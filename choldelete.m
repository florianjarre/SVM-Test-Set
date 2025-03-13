function Rt = choldelete(R,d)
% Given a nonsingular upper triangular matrix R with R'*R = A, 
% let At be obtained from A by deleting row d and column d. 
% Then compute Rt with Rt'*Rt = At.

[~,n] = size(R);
if d < 1 || d > n
    error('dimensions in choldelete do not fit');
end
Rt= zeros(n-1);
if d>1  % this if-query is not really necessary....
    Rt(1:d-1,1:d-1) = R(1:d-1,1:d-1);
    Rt(1:d-1,d:n-1) = R(1:d-1,d+1:n);
end
Rt(d:n-1,d:n-1) = cholupdate(R(d+1:n,d+1:n),R(d,d+1:n).');
