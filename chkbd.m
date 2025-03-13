function[X,z,x1low,x1high,x2low,x2high] = chkbd(nnn,m)
% generate SVM data in 2d-checker board pattern
% Generate a matrix with columns x^i (1\le i\le m)
% and a vector z with components \pm 1 classifying the vectors x^i
% Also specify rangen in which a plot is to be made
%
% Input: for example:
% nnn = 3;  % checker board of size nnn times nnn 
% m = 500;  % m points on a nnn by nnn checker board

X = nnn*rand(2,m);

Xlow = [floor(X(1,:));floor(X(2,:))];
z = sign(mod(mod(Xlow(1,:),2)+mod(Xlow(2,:),2) ,2)-.5).';

% indicate the domain, where the plot is to be done
% x1low = -1; x1high = 4;
% x2low = -1; x2high = 4;
x1low = 0; x1high = 3; % for plotting
x2low = 0; x2high = 3; % for plotting

end