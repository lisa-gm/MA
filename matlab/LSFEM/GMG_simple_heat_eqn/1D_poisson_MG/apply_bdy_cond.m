% ******************************************************************* %
% ***************** APPLYING BDY CONDITIONS TO RHS ****************** %
% ******************************************************************* %

function [rhs_sym, Hess_J_sym] = apply_bdy_cond(Nx_pts, Nt_pts, Hess_J, f_vec, bdy_cond, u0, bdy_left, bdy_right, level)

if(level == 1)
tot_pts = Nx_pts*Nt_pts;

rhs_sigma = f_vec(1:tot_pts);
rhs_u = f_vec(tot_pts+1:end);

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

bdy_ind = get_bdy_ind(Nx_pts, Nt_pts, bdy_cond);

Hess_J_sym = Hess_J;
Hess_J_sym(bdy_ind,:) = 0;

%% symmetrize the matrix

rhs_sym(:) = rhs(:)- Hess_J_sym(:,bdy_ind)*rhs(bdy_ind);
rhs_sym = rhs_sym';

Hess_J_sym(:,bdy_ind) = 0;

%% put in ones in diagonal

tempA = speye(length(Hess_J));
tempA(bdy_ind,:) = 0;
tempB = speye(length(Hess_J));
tempD = tempB-tempA;

Hess_J_sym = Hess_J_sym + tempD;

else
tot_pts = Nx_pts*Nt_pts;

bdy_ind = get_bdy_ind(Nx_pts, Nt_pts, bdy_cond);

rhs_sym = zeros(2*tot_pts,1);
% rhs doesnt matter here ... 

Hess_J_sym = Hess_J;
Hess_J_sym(bdy_ind,:) = 0;
Hess_J_sym(:,bdy_ind) = 0;

%% put in ones in diagonal

tempA = speye(length(Hess_J));
tempA(bdy_ind,:) = 0;
tempB = speye(length(Hess_J));
tempD = tempB-tempA;

Hess_J_sym = Hess_J_sym + tempD;
    
end
end