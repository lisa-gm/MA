% *********************************************************************** %
% ********************** BLOCK JACOBI SMOOTHER ************************** %
% *********************************************************************** %

% considering 2x2 blocks of sigma and u for same dof

function [u, sol] = JacobiSolve_LS(A, f, u0, max_iter)
tot_pts = length(f)/2;
% additional values, convergence threshold eps,
% damping factor hard coded for now
omega = 0.6656;
%omega = 1;
eps = 10^(-12);

u_old = u0;
u_new = zeros(size(u0));

sol = zeros(length(f), max_iter);

for iter = 1:max_iter
    for i=1:tot_pts
        c = A([i,i+tot_pts],[i,i+tot_pts]) \ (f([i, i+tot_pts]) - A([i,i+tot_pts], :)*u_old);
        %u_old
        %rhs = (f([i, i+tot_pts]) - A([i,i+tot_pts], :)*u_old)
      
        u_new([i, i+tot_pts]) = u_old([i, i+tot_pts]) + omega*c;
    end
    
    u_old = u_new;
    sol(:, iter) = u_new;
    if(norm(f-A*u_new) < eps)
       sol = sol(:, 1:iter);
       break;
    end
    
    %fprintf('norm res : %d\n', norm(f-A*u_new));
    
end

u = u_new;
end