%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% FUNCTION TO DO A DIRECT SOLVE %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [sigma_direct, u_direct] = direct_solve_lsfem(Nx_pts, Nt_pts, inner_ind_sigma, inner_ind_u, bdy_cond, bdy_left, bdy_right, u0, Hess_J, grad_J, sigma, u) 
tot_pts = Nx_pts*Nt_pts;
update_direct = zeros(2*tot_pts,1);

update_direct = - Hess_J \ grad_J;
%update_direct([inner_ind_sigma, tot_pts+inner_ind_u])=update_inner_direct;

update_direct_sigma = update_direct(1:tot_pts);
update_direct_u = update_direct(tot_pts+1:end);

sigma_direct = sigma + update_direct_sigma;
u_direct = u + update_direct_u;

u_mat_direct = zeros(Nx_pts, Nt_pts);
sigma_mat_direct = zeros(Nx_pts, Nt_pts);

for ts=1:Nt_pts
    u_mat_direct(:, ts) = u_direct((ts-1)*Nx_pts+1:ts*Nx_pts);
    sigma_mat_direct(:, ts) = sigma_direct((ts-1)*Nx_pts+1:ts*Nx_pts);
end

if(strcmp(bdy_cond, 'Dirichlet'))
    u_mat_direct(:,1) = u0;
    u_mat_direct(1, :) = bdy_left;
    u_mat_direct(end, :) = bdy_right;
end

if(strcmp(bdy_cond, 'Neumann'))
    u_mat_direct(:,1) = u0;
    sigma_mat_direct(1, :) = bdy_left;
    sigma_mat_direct(end, :) = bdy_right;
end

sigma_direct = reshape(sigma_mat_direct, [Nx_pts*Nt_pts,1]);
u_direct = reshape(u_mat_direct, [Nx_pts*Nt_pts,1]);

end