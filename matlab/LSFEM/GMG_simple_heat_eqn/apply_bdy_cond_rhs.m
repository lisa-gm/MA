% ******************************************************************* %
% ***************** APPLYING BDY CONDITIONS TO RHS ****************** %
% ******************************************************************* %

function rhs = apply_bdy_cond(H, f_sigma_vec, f_u_vec, bdy_cond, u0, bdy_left, bdy_right)

rhs_sigma = f_sigma_vec;
rhs_u = f_u_vec;

rhs_u(1:Nx_pts) = u0;

if(strcmp(bdy_cond, 'Dirichlet'))
    rhs_u(1:Nx_pts:end) = bdy_left;
    rhs_u(Nx_pts:Nx_pts:end) = bdy_right;
end

if(strcmp(bdy_cond, 'Neumann'))
    rhs_sigma(1:Nx_pts:end) = bdy_left;
    rhs_sigma(Nx_pts:Nx_pts:end) = bdy_right;
end

rhs = [rhs_sigma; rhs_u];

Hess_J_sym = Hess_J;
Hess_J_sym(bdy_ind,:) = 0;

%% symmetrize the matrix

rhs_sym(:) = rhs(:) - Hess_J_sym(:,bdy_ind)*rhs(bdy_ind);
rhs_sym = rhs_sym';
Hess_J_sym(:,bdy_ind) = 0;


end