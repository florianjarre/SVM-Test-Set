% kkt_norm
% compute some relative norm of the violation of the KKT conditions 
% for the problem
%    min 0.5*x'*H*x -c'*x | z'*x = 0; C*e >= x >= 0
% at a given point x.

tolAS = 1.0e-14; % ignore error of the (relative) norm tolAS when
                 % determining the active set.

kviol1 = -min(x);
kviol1 = max(0,kviol1);
kviol2 = max(x-C);
kviol2 = max(0,kviol2);
kviol3 = abs(z.'*x)/(norm(x)+eps); % "+eps" to allow x = 0 
ACT1 = x <= tolAS;
ACT2 = x >= (1-tolAS)*C;
INACT = ~(ACT1|ACT2);
nIN = sum(INACT);
g = H*x-c;
if nIN > 0
    mu = g(INACT).'*z(INACT)/nIN;
    gt = g-mu*z; % KKT implies gt(ACT1) >= 0,  gt(ACT2) <= 0
    if sum(ACT1) > 0
        gt(ACT1) = min(0,gt(ACT1));
    end
    if sum(ACT2) > 0
        gt(ACT2) = max(0,gt(ACT2));
    end
    kviol4 = norm(gt)/norm(x); % relative norm of projected gradient
else
    sigma = zeros(size(x));
    sigma(ACT1) = 1;
    sigma(ACT2) = -1;
    tmp1 = sigma.*z > 0;
    tmp2 = sigma.*z < 0;
    kviol4 = min(sigma(tmp1).*g(tmp1))-max(-sigma(tmp2).*g(tmp2));
    kviol4 = max(0,-kviol4);
end
violkkt = max([kviol1,kviol2,kviol3,abs(kviol4)]);

