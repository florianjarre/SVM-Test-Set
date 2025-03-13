% test the halfmoon problem with sweep, gsmo, and rsmo

d = 3;delta = .25;rho = 0.25; n = 250;
% d:     dimension of data space,         (should be less than 50)
% delta: width of half moon,              (should be about 0.25)
% rho:   shift of S_2,                    (should be >= delta)
% n:     number of training data points   (should be less than 10000)

[X,z] = halfmoon(d,delta,rho,n);
[H,C,gamma] = gen_kernel(X,z);
c = ones(size(z));

disp('run sweep')
tic;[x,outfile] = sweep(H,z,C);outfile,toc
disp('relative classification errors - false negative, false positive')
testhalfmoon

disp('now run gsmo')
tic;[x,outfile] = gsmo(H,z,C);outfile,toc
disp('relative classification errors - false negative, false positive')
testhalfmoon

disp('now run rsmo')
tic;[x,outfile] = rsmo(H,z,C);outfile,toc
disp('relative classification errors - false negative, false positive')
testhalfmoon
