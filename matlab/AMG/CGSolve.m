%%% CONJUGATE GRADIENT METHOD %%%%%%

function [u, x, res, r] = CGSolve(A, f, x0, n)
eps = 10^(-10);

dim = length(x0);
p = zeros(dim,n+1);
r = zeros(dim,n+1);
x = zeros(dim,n+1);

p(:,1) = f - A*x0;
r(:,1) = p(:,1);
x(:,1) = x0;
res = zeros(1, n+1);
res(1) = norm(r(:,1));

% saving all r,x,p values for statistics later
    for i=2:n+1
    a = (r(:,i-1)'*r(:,i-1))/dot(p(:,i-1),A*p(:,i-1));
    x(:,i) = x(:,i-1) + a*p(:,i-1);
    r(:,i) = r(:,i-1) - a*A*p(:,i-1);
    
    if(norm(r(:,i)) < eps)
        x = x(:,1:i);
        res(i) = norm(r(:,i));
        res = res(1:i);
        r = r(:,1:i);
        % fprintf('converged after %d steps', i);
        break; 
    end
    
    b = (r(:,i)'*r(:,i))/(r(:,i-1)'*r(:,i-1));
    p(:,i) = r(:,i) + b*p(:,i-1);
    res(i) = norm(r(:,i));
    end
    
% only solution vector
u = x(:,end);
end