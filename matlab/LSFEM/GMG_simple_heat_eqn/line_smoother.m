% ******************************************************************** %
% ***************** LINE SMOOTHER -- 2 COLOURING ********************* %
% ******************************************************************** %

% extract all entries belonging to spatial position x_i, direct solve
% for all time steps and \sigma and u

% but go in the order x_1, x_3, x_5, ... and then
% update solution, recompute residual and then solve for
% x_2, x_4, ... 

% dont do anything special on bdy, should be in the rhs
%  

function [sol, all_sol] = line_smoother(Nx_pts, Nt_pts, A, rhs, sol, max_iter, eps)
tot_pts = Nx_pts*Nt_pts;

% determine sets of indices
% each column belongs to indiceus of one time steps
ind_sp_odd = transpose(1:2:Nx_pts);
len_sp_odd = length(ind_sp_odd);

ind_sp_even = transpose(2:2:Nx_pts);
len_sp_even = length(ind_sp_even);

ind_odd = zeros(len_sp_odd, Nt_pts);
ind_even = zeros(len_sp_even, Nt_pts);

all_sol = zeros(length(rhs), max_iter);

for ts = 1:Nt_pts
    ind_odd(:, ts) =  ind_sp_odd + (ts-1)*Nx_pts;
    ind_even(:, ts) =  ind_sp_even + (ts-1)*Nx_pts;   
end

% for sigma & u
ind_odd = horzcat(ind_odd, ind_odd+tot_pts);
ind_even = horzcat(ind_even, ind_even+tot_pts);

for it = 1:max_iter
    
    res = rhs - A*sol;
    
    % extract & solve for "lines", ODD
    for ind=1:length(ind_sp_odd)
        subM = A(ind_odd(ind,:), ind_odd(ind,:));
        subRes = res(ind_odd(ind,:));

        sub_update = subM \ subRes;
        sol(ind_odd(ind,:)) = sol(ind_odd(ind, :)) + sub_update;
    end

    res = rhs - A*sol;
    % extract & solve for "lines", EVEN
    for ind=1:length(ind_sp_even)
        subM = A(ind_even(ind,:), ind_even(ind,:));
        subRes = res(ind_even(ind,:));

        sub_update = subM \ subRes;
        sol(ind_even(ind,:)) = sol(ind_even(ind, :)) + sub_update;
    end
    
     all_sol(:, it) = sol;
     
     if(norm(rhs-A*sol) < eps)
      fprintf('line smoother converged with res %d\n', norm(A*sol - rhs));
       all_sol = all_sol(:, 1:it);
       break;
     end
    
        u_ls = sol(tot_pts+1:end);
        u_mat_ls = zeros(Nx_pts, Nt_pts);

        for ts=1:Nt_pts
            u_mat_ls(:, ts) = u_ls((ts-1)*Nx_pts+1:ts*Nx_pts);
        end
        
    % fprintf('res after %d iter: %d\n', it, norm(A*sol-rhs));
end


end