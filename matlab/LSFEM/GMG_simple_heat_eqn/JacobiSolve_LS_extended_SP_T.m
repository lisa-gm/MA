% *********************************************************************** %
% ***************** BLOCK JACOBI SMOOTHER EXTENDED! ********************* %
% *********************************************************************** %

% considering blocks of size block_size for sigma and u for same dof on
% grid

% here we walk in time then space, could also couple time or do patches 
% of space-time

function [u, sol] = JacobiSolve_LS_extended_SP_T(Nx_pts, Nt_pts, A, f, u0, max_iter, omega, block_size_s, block_size_t)
tot_pts = length(f)/2;

eps = 10^(-9);

u_old = u0;
u_new = zeros(size(u0));

sol = zeros(length(f), max_iter);

% since block size can vary maybe not divisible
% careful to not couple grid(i, end) and grid(i+1, 1) ... 
k_s = floor(Nx_pts/block_size_s);
last_mult_s = k_s*block_size_s;
rest_size_s = Nx_pts - last_mult_s;

k_t = floor(Nt_pts/block_size_t);
last_mult_t = k_t*block_size_t;
rest_size_t = Nt_pts - last_mult_t;


for iter = 1:max_iter
    for j = 1:block_size_t:last_mult_t
        for i=1:block_size_s:last_mult_s
            % get list of sigma and u indices we want
            %patch_ind = [i+(j-1)*Nx_pts,i+(j-1)*Nx_pts+tot_pts];
            %c_small = A(patch_ind, patch_ind) \ (f(patch_ind) - A(patch_ind, :)*u_old)
            
            curr_inds = [];
            for l=0:block_size_t-1
                curr_inds = [curr_inds, (j-1+l)*Nx_pts+i: (j-1+l)*Nx_pts+i-1+block_size_s, (j-1+l)*Nx_pts+i+tot_pts: (j-1+l)*Nx_pts+i-1+block_size_s+tot_pts];   
            end
            
            
            curr_inds = sort(curr_inds);
            
            c = A(curr_inds, curr_inds)  \ (f(curr_inds) - A(curr_inds, :)*u_old);
   
            
            u_new(curr_inds) = u_old(curr_inds) + omega*c;
        end
        
        % rest size indices never reached ... for specific t, update now
        if(rest_size_s ~= 0)
            curr_inds = [];
            
            for l=0:block_size_t-1
                curr_inds = [curr_inds, (j-1+l)*Nx_pts+[last_mult_s+1:Nx_pts], tot_pts+(j-1+l)*Nx_pts+[last_mult_s+1:Nx_pts]];
            end
            curr_inds = sort(curr_inds);
            
            c = A(curr_inds, curr_inds)  \ (f(curr_inds) - A(curr_inds, :)*u_old);
            u_new(curr_inds) = u_old(curr_inds) + omega*c;
        end
          
    end
    
    if(rest_size_t ~= 0)

        
    for i=1:block_size_s:last_mult_s
    curr_inds = [];
        for l=0:rest_size_t-1
            curr_inds = [curr_inds, (last_mult_t+l)*Nx_pts+i: (last_mult_t+l)*Nx_pts+i-1+block_size_s, (last_mult_t+l)*Nx_pts+i+tot_pts: (last_mult_t+l)*Nx_pts+i-1+block_size_s+tot_pts];
        end
    curr_inds = sort(curr_inds);
            
    c = A(curr_inds, curr_inds)  \ (f(curr_inds) - A(curr_inds, :)*u_old);
    u_new(curr_inds) = u_old(curr_inds) + omega*c;
    end
        
        
    curr_inds = [];
            
    for l=0:rest_size_t-1
        curr_inds = [curr_inds, (last_mult_t+l)*Nx_pts+[last_mult_s+1:Nx_pts], tot_pts+(last_mult_t+l)*Nx_pts+[last_mult_s+1:Nx_pts]];
    end
    curr_inds = sort(curr_inds);
            
    c = A(curr_inds, curr_inds)  \ (f(curr_inds) - A(curr_inds, :)*u_old);
    u_new(curr_inds) = u_old(curr_inds) + omega*c;
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