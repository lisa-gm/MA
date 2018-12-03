function [u, error] = JacobiSolve(A, f, u0, max_iter)

% additional values, convergence threshold eps,
% damping factor
eps = 10^(-12);
w = 1/2;
N = length(f);

% get Jacobi matrix
invJ = diag(1./(diag(A)));
u_new = zeros(size(u0));

% first iterate
% use additional damping factor, can be just 1
% u_new = u0 + w*invJ*(f- A*u0);
  for i=1:N   
    sum_1 = sum(A(i,1:end).*u0');
    u_new(i) = u0(i) + w*invJ(i)*(f(i)- sum_1);
  end

u_new_1 = u0 + w*invJ*(f - A*u0);
% save solutions u_i
error = zeros(1,max_iter+1);
error(1) = norm(f- A*u0);
error(2) = norm(f- A*u_new);

for k=2:max_iter
    if(norm(error(k)) < eps)
        %fprintf('Jacobi iteration converged, iter count: %d\n', i);
        break;
    end
    
    % elementwise implementation so that we dont need inverse
    for i=1:N   
    sum_1 = sum(A(i,1:end).*u0');
    u_new(i) = u0(i) + w*invJ(i)*(f(i)- sum_1);
    end
    u_new = u0 + w*invJ*(f - A*u0);
    
    u0 = u_new;
    error(k+1) = norm(f- A*u_new);
end

u = u_new;
end