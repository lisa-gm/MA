% *********************************************************************** %
% ***************** BLOCK JACOBI SMOOTHER EXTENDED! ********************* %
% *********************************************************************** %

% considering blocks of size block_size for sigma and u for same dof on
% grid

% here we walk in time then space, could also couple time or do patches 
% of space-time

function [u, sol] = JacobiSolve_LS_extended(Nx_pts, Nt_pts, A, f, u0, max_iter, omega, block_size)
tot_pts = length(f)/2;

% additional values, convergence threshold eps,
% damping factor hard coded for now
%omega = 0.6656;
%omega = 1;
eps = 10^(-12);

u_old = u0;
u_new = zeros(size(u0));

sol = zeros(length(f), max_iter);

% since block size can vary maybe not divisible
% careful to not couple grid(i, end) and grid(i+1, 1) ... 
k = floor(Nx_pts/block_size);
last_mult = k_s*block_size;
rest_size = Nx_pts - last_mult;

for iter = 1:max_iter
    for j = 1:Nt_pts
        for i=1:block_size_s:last_mult_s
            % get list of sigma and u indices we want
            %patch_ind = [i+(j-1)*Nx_pts,i+(j-1)*Nx_pts+tot_pts];
            %c_small = A(patch_ind, patch_ind) \ (f(patch_ind) - A(patch_ind, :)*u_old)
            
            for l=1:block_size_t
                curr_inds = [(j-1)*Nx_pts+i: (j-1)*Nx_pts+i-1+block_size, (j-1)*Nx_pts+i+tot_pts: (j-1)*Nx_pts+i-1+block_size+tot_pts];
                
            end
            c = A(curr_inds, curr_inds)  \ (f(curr_inds) - A(curr_inds, :)*u_old);
            %u_old
            %rhs = (f(curr_inds) - A(curr_inds, :)*u_old)
            
            u_new(curr_inds) = u_old(curr_inds) + omega*c;
        end
        
        % rest size indices never reached ... for specific t, update now
        if(rest_size ~= 0)
            curr_inds = [(j-1)*Nx_pts+[last_mult+1:Nx_pts], tot_pts+(j-1)*Nx_pts+[last_mult+1:Nx_pts]];
            c = A(curr_inds, curr_inds)  \ (f(curr_inds) - A(curr_inds, :)*u_old);
            u_new(curr_inds) = u_old(curr_inds) + omega*c;
        end
       
            
    
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