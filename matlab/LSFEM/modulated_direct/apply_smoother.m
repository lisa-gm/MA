% *********************************************************************** %
% ***************** APPLY GIVEN SMOOTHER! ********************* %
% *********************************************************************** %

 

function [u, sol] = apply_smoother(A, P, f, u0, omega, max_iter)
eps = 10^(-9);

sol = zeros(length(u0), max_iter+1);
sol(:,1) = u0;

for it=1:max_iter
    
    u0 = u0 + omega*P*(f - A*u0);
    sol(:, it+1) = u0;  
    
    if(norm(f - A*u0) < eps)
        sol = sol(:,it+1);
        fprintf('block Jacobi smoother converged after %d iter', it);
        break;
    end
        
end

u = u0;    
end