function [A,M]=AssembleMA(n)
h = 1/(n-1);

e = h*ones(n,1);
M = spdiags([e/6 2*e/3 e/6],-1:1,n,n);
M(1,1) = h/3;
M(n,n) = h/3;

e = (1/h)*ones(n,1);
A = spdiags([-e 2*e -e],-1:1,n,n);
A(1,1) = 1/h;
A(n,n) = 1/h;
return