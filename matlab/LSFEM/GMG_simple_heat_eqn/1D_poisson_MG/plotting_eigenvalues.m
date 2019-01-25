% ********************************************************************* %
% ********************* plot eigenvalues to compare ******************** %
% ********************************************************************* %
clear all;
close all;

% CHOOSE GENERAL PARAMTERS: 

% MG parameters
max_iter = 200;                          % number of newton iterations
eps = 10^(-9);                           % tolerance
levels = 2; 
smoother = 'Jacobi_LS';                  

% eqn parameters
diff_const = 1;                          % diffusion constant, define through sigma = diff_const*nabla(u) 

% LSFEM parameters
c1 = 2;                                
c2 = 2;                               
            
% ********** SETTING UP DOMAIN PARAMETERS ***************************** %

% been using 12 x 12 points in space time so 144 points here  

% space 
S = 1;
Nx_elem = 144;
Nx_pts = Nx_elem + 1;

x_vec_f = linspace(0,S, Nx_pts);

% ******** BOUNDARY CONDITIONS ON FINEST GRID ************************** %

bdy_left = 1;
bdy_right = 0;

% ***************** INITIALISE U & SIGMA w/ bdy ******************************* %
sigma = zeros(Nx_pts,1);
u = zeros(Nx_pts,1); 
u(1) = bdy_left;
u(2) = bdy_right;

sol_init = [sigma; u];

% ********************* CONSTRUCTING MATRICES *********************** %

% same as reaction diffusion, just took out u_t part
[J_ss_lin, J_su_lin, J_us_lin, J_uu_lin] = poisson_mat_no_bdy_lsfem (S, Nx_elem, diff_const);
Hess_J = [J_ss_lin, J_su_lin; J_us_lin, J_uu_lin];


% ****************** RHS ******************************************** %

f_sigma_vec = zeros(Nx_pts,1);
f_u_vec = zeros(Nx_pts,1);
%f_u_vec = transpose(x_vec_f);
f = [f_sigma_vec; f_u_vec];

rhs = poisson_rhs_no_bdy_lsfem (S, Nx_elem, f);

rhs(Nx_pts+1) = bdy_left;
rhs(end) = bdy_right;

[rhs_sym, Hess_J_sym] = apply_bdy_cond_poisson(Nx_pts, Hess_J, rhs, bdy_left, bdy_right, 1);

% ********************** JACOBI ITER ********************************** %

[sol_Jac, sol] = JacobiSolve_LS(Hess_J_sym, rhs_sym, sol_init, 300);

sigma_Jac = sol_Jac(1:Nx_pts);
u_Jac = sol_Jac(Nx_pts+1:end);

% ************* DIRECT SOLVE FOR COMPARISON ************************* %

sol_direct = Hess_J_sym \ rhs_sym;

sigma_direct = sol_direct(1:Nx_pts);
u_direct = sol_direct(Nx_pts+1:end);

% ******************** PLOTTING EIGENVALUES ******************** %

[V_sym, D_sym] = eig(Hess_J_sym);
D_sym = sort(diag(D_sym));

figure;
scatter(1:length(D_sym), D_sym, '.');
ylabel('eigenvalues', 'interpreter','latex');

fprintf('condition number symmetrised hessian %d\n', cond(Hess_J_sym));

% figure;
% plot(x_vec_f, sol_direct(Nx_pts+1:end));
