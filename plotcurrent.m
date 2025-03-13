% plotcurrent.m
% after solving the dual kernel problem determine the associated value
% beta and then plot the regions divided by the svm
% Input: plotting range   x1low,x1high,x2low,x2high
% and                     H, c, z and the solution x from
%    min 0.5*x'*H*x -c'*x | z'*x = 0; C*e >= x >= 0.
% ( C is not used )
 
nplot = 200; % number of points (minus one) in x1- and x2-direction 

n = length(z); % adjust notations....
c = ones(n,1);

xplot1 = 0:nplot;
xplot2 = xplot1;
onp = ones(nplot+1,1)/nplot;
XR1 = x1low+(x1high-x1low)*xplot1.'*onp.';
XR2 = x2low+(x2high-x2low)*onp*xplot2;
XR = [XR1(:).';XR2(:).'];
spp = length(XR1(:)); % number of sample points
zR = zeros(1,spp);

beta = findbeta(x,H,z,C);

for ii = 1:spp
xtilde = XR(:,ii);
kappatx = 0;
for i = 1:n
    kappatx = kappatx + c(i)*x(i)*z(i)*exp(-gamma*norm(X(:,i)-xtilde)^2);
end
ztilde = sign(kappatx-beta);
zR(ii) =  ztilde;
end

I1 = zR > 0;
XR1 = XR(:,I1);
plot(XR1(1,:),XR1(2,:),'g.')    % green area
hold on
I2 = zR < 0;
XR2 = XR(:,I2);
plot(XR2(1,:),XR2(2,:),'r.')    % red area
plot(X(1,:),X(2,:),'k*')        % data points 

% plot the grid next
for i = ceil(x1low):floor(x1high)
    plot([x1low,x1high],[i,i],'b-')
end
for j = ceil(x2low):floor(x2high)
   plot([j,j],[x2low,x2high],'b-')
end

hold off
ax = gca;ax.FontSize = 15; % make numbers at the axis bigger
