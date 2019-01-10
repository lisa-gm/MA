% ******************************************************************* %
% ***************** APPLYING BDY CONDITIONS TO RHS ****************** %
% ******************************************************************* %

function [rhs_sym, Hess_J_sym] = apply_bdy_cond_poisson(Nx_pts, Hess_J, f_vec, bdy_left, bdy_right, level)

if(level == 1)
rhs_sigma = f_vec(1:Nx_pts);
rhs_u = f_vec(Nx_pts+1:end);

rhs_u(1) = bdy_left;
rhs_u(end) = bdy_right;

rhs = [rhs_sigma; rhs_u];

bdy_ind = [Nx_pts+1, 2*Nx_pts];

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
bdy_ind = [Nx_pts+1, 2*Nx_pts];

rhs_sym = zeros(2*Nx_pts,1);
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