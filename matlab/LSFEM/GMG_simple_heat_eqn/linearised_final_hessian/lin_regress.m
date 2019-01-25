% ******************************************************** %
% ********* lineares auslgeichsproblem loesen ************ % 
% ******************************************************** %

% include k=0 ?!
function [c, rho] = lin_regress(no_iter, norm_res)

% construct matrix A
A = zeros(no_iter+1, 2);
A(:,1) = 1;
A(:,2) = 0:no_iter;

b = transpose(log(norm_res));

% x = [log(c); log(rho)]
x = (A'*A) \ (A'*b);

c = exp(x(1));
rho = exp(x(2));
end