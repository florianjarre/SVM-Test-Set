Matlab Test Set for Support Vector Machines with Gaussian Kernel


Generating test data,  

a matrix X whose colums are the data points and a column vector z with associated labels
                       
via

   [X,z,x1low,x1high,x2low,x2high] = chkbd(nnn,m);
   
or

   [X,z] = halfmoon(d,delta,rho,n);
   
(In chkbd.m also the range for the plot is given as output

and can be used in plotcurrent.m after solving the SVM problem)


Making the kernel via

   [H,C,gamma] = gen_kernel(X,z);
   
(with constants specified in gen_kernel.m)


Solving the SVM problem via

   [x,outfile] = sweep(H,z,C);
   
or

   [x,outfile] = gsmo(H,z,C);
   
or

   [x,outfile] = rsmo(H,z,C);
   

Sub-programs used are

kkt_norm.m

pilar.m

findbeta.m


testhalfmoon.m tests the computed classification of the halfmoon set

plotcurrent.m plots the checkerborad classification
