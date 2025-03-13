function J = pilar(w,I,m)
% PIck m LARgest elements of w(I)
if m > sum(I)
    error('too few elements to pick from in pilar')
end
n = length(w);
w(~I) = -Inf;
[~,tmp] = sort(w,'descend');
J = false(n,1);
J(tmp(1:m)) = true;
end

