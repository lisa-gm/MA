% *********************************************************************** %
% ********************** BLOCK JACOBI SMOOTHER ************************** %
% *********************************************************************** %

% considering 2x2 blocks of sigma and u for same dof

function [u, res_list] = JacobiSolve_LS(A, f, u0, max_iter)
tot_pts = length(f)/2;
% additional values, convergence threshold eps,
% damping factor hard coded for now
omega = 0.6656;
%omega = 1;
eps = 10^(-12);

u_old = u0;
u_new = zeros(size(u0));

res_list = zeros(1, max_iter+1);
res_list(1) = norm(f-A*u_old);

for iter = 1:max_iter
    for i=1:tot_pts
        c = A([i,i+tot_pts],[i,i+tot_pts]) \ (f([i, i+tot_pts]) - A([i,i+tot_pts], :)*u_old);
        %u_old
        %rhs = (f([i, i+tot_pts]) - A([i,i+tot_pts], :)*u_old)
      
        u_new([i, i+tot_pts]) = u_old([i, i+tot_pts]) + omega*c;
    end
    
    u_old = u_new;
    res = norm(f-A*u_new);
    res_list(1, iter+1) = res;
    if(res < eps)
       fprintf('Jacobi_LS converged after %d iter,with res %d\n', iter, res); 
       res_list = res_list(1, 1:iter+1);
       break;
    end
    
    %fprintf('norm res : %d\n', norm(f-A*u_new));
    
end

u = u_new;
end