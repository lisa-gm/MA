% ******************************************************** %
% *********** MULTIGRID AROUND THE SOLUTION ************** % 
% ******************************************************** %

clear all;
close all;

% SOLVE H*u = r for r = H*u_old - grad_J;

% load hessian
% contains the variables 'Hess_J', 'grad_J', 'S', 'Nx_elem', 'T', 'Nt_elem','sigma', 'u', 'u0', 'bdy_cond', 'bdy_left', 'bdy_right','diff_const')
load('Hess_J_20by40elem_diff0_001.mat');

% MG parameters
max_iter = 300;                          % number of newton iterations
eps = 10^(-9);                           % tolerance
levels = 2; 
% smoother = 'line_smoother';
smoother = 'Jacobi_LS_extended_SP_T';

% eqn parameters
diff_const = 1;                       

% LSFEM parameters
c1 = 2;                               
c2 = 2; 

% OTHER PARAMETERS
Nx_pts = Nx_elem + 1;
Nt_pts = Nt_elem + 1;
tot_pts = Nx_pts*Nt_pts;

x_vec = linspace(0,S, Nx_pts);
t_vec = linspace(0,T, Nt_pts);

Nx_elem_c = Nx_elem/2;
Nt_elem_c = Nt_elem/2;

Nx_elem_list = [Nx_elem, Nx_elem_c];
Nt_elem_list = [Nt_elem, Nt_elem_c];

Nx_pts_list = Nx_elem_list + 1;
Nt_pts_list = Nt_elem_list + 1;

hx_list = [S/Nx_elem, S/Nx_elem_c];
ht_list = [T/Nt_elem, T/Nt_elem_c]; 

% ***************** INITIALISE U & SIGMA w/ bdy ******************************* %

sigma_bdy_mat = ones(Nx_pts_list(1), Nt_pts_list(1));
u_bdy_mat = zeros(Nx_pts_list(1), Nt_pts_list(1));

for ts=1:Nt_pts_list(1)
    u_bdy_mat(:,ts) = u0;
end

if(strcmp(bdy_cond, 'Dirichlet'))
    u_bdy_mat(1, :) = bdy_left;
    u_bdy_mat(end, :) = bdy_right;
end 

if(strcmp(bdy_cond, 'Neumann'))
    sigma_bdy_mat(1, :) = bdy_left;
    sigma_bdy_mat(end, :) = bdy_right;
end

sigma_bdy = reshape(sigma_bdy_mat, [tot_pts,1]);
u_bdy = reshape(u_bdy_mat, [tot_pts,1]);

sol_init = [sigma_bdy; u_bdy];

% construct rhs here
r = Hess_J*[sigma; u] - grad_J;

% apply boundary conditions
[rhs_sym, Hess_J_sym] = apply_bdy_cond(Nx_pts, Nt_pts, Hess_J, r, bdy_cond, u0, bdy_left, bdy_right, 1);

% ********************** MULTIGRID ********************************** %

[sol_mg, it_sol_mg, it_res_mg] = V_cycle_galerkin_assembly(Nx_elem_list, hx_list, Nt_elem_list, ht_list, Hess_J_sym, rhs_sym, sol_init, bdy_cond, u0, bdy_left, bdy_right, levels, max_iter, smoother, eps);

sigma_mg = sol_mg(1:tot_pts);
u_mg = sol_mg(tot_pts+1:end);

% ************* DIRECT SOLVE FOR COMPARISON ************************* %

sol_direct = Hess_J_sym \ rhs_sym;

sigma_direct = sol_direct(1:tot_pts);
u_direct = sol_direct(tot_pts+1:end);


% ************************ PLOTTING ******************************** %

% reshape
u_mat_direct = zeros(Nx_pts_list(1), Nt_pts_list(1));
u_mat_mg = zeros(Nx_pts_list(1), Nt_pts_list(1));

sigma_mat_direct = zeros(Nx_pts_list(1), Nt_pts_list(1));
sigma_mat_mg = zeros(Nx_pts_list(1), Nt_pts_list(1));


for ts=1:Nt_pts_list(1)
    sigma_mat_direct(:,ts) = sigma_direct((ts-1)*Nx_pts_list(1)+1:ts*Nx_pts_list(1));
    u_mat_direct(:, ts) = u_direct((ts-1)*Nx_pts_list(1)+1:ts*Nx_pts_list(1));
    
    %sigma_mat_direct(:,ts) = sigma((ts-1)*Nx_pts_list(1)+1:ts*Nx_pts_list(1));
    %u_mat_direct(:, ts) = u((ts-1)*Nx_pts_list(1)+1:ts*Nx_pts_list(1));    
    
    sigma_mat_mg(:, ts) = sigma_mg((ts-1)*Nx_pts_list(1)+1:ts*Nx_pts_list(1));
    u_mat_mg(:, ts) = u_mg((ts-1)*Nx_pts_list(1)+1:ts*Nx_pts_list(1));
end

figure;
mesh(x_vec, t_vec, transpose(u_mat_direct)); 
xlabel('space');
ylabel('time');
zlabel('u');
%zlim([min(u), max(u)]);
title('u direct');

figure;
mesh(x_vec, t_vec, transpose(u_mat_mg)); 
xlabel('space');
ylabel('time');
zlabel('u');
%zlim([min(u), max(u)]);
title('u multigrid');