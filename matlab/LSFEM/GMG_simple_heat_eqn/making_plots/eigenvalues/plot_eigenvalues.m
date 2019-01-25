% ******************************************************************* %
% ********* CHECKING FOR EIGENVALUES & DIAGONAL DOMINANCE *********** %
% ******************************************************************* %


% construct and then save Hessian for 11 x 11 
% and 24 x 44, check eigenvalues etc. 
% potentially also for iteration matrix of block Jacobi 1x1

clear all;
close all;

% BLOCK SIZE
block_size_s = 1;
block_size_t = 1;

% eqn parameters
diff_const = 0.1;                          % diffusion constant, define through sigma = diff_const*nabla(u) 

% LSFEM parameters
c1 = 2;                                 % constants in front of J = c1*|| ....-f(u) ||^2 + c2*|| ...||^2
c2 = 2;                                 % need to multiply with them in computation of each integral 
                                        
% ********** SETTING UP DOMAIN PARAMETERS ***************************** %
% space 
S = 1;
Nx_elem = 24;
Nx_pts = Nx_elem + 1;

% time
T = 1;
Nt_elem = 44;
Nt_pts = Nt_elem + 1;

tot_pts = Nx_pts*Nt_pts;

% ************************ BOUNDARY CONDITIONS ************************** %

%bdy_cond = 'Dirichlet';
bdy_cond = 'Neumann';

u0 = zeros(Nx_pts,1);

bdy_left = zeros(Nt_pts,1);
bdy_right = zeros(Nt_pts,1);


% ***************** INITIALISE U & SIGMA w/ bdy ******************************* %

% create vector with random variables between [0,1] as initial guess
sigma_init_mat = rand(Nx_pts, Nt_pts);
u_init_mat = rand(Nx_pts, Nt_pts);

u_init_mat(:,1) = u0;

if(strcmp(bdy_cond, 'Dirichlet'))
    u_init_mat(1, :) = bdy_left;
    u_init_mat(end, :) = bdy_right;
end 

if(strcmp(bdy_cond, 'Neumann'))
    sigma_init_mat(1, :) = bdy_left;
    sigma_init_mat(end, :) = bdy_right;
end

sigma_init = reshape(sigma_init_mat, [tot_pts,1]);
u_init = reshape(u_init_mat, [tot_pts,1]);

sol_init = [sigma_init; u_init];

% ****************** CONSTRUCTING MATRICES ****************** %

[J_ss_lin, J_su_lin, J_us_lin, J_uu_lin] = lin_part_mat_no_bdy_lsfem (S, Nx_elem, T, Nt_elem, c1 ,c2, diff_const);
Hess_J = [J_ss_lin, J_su_lin; J_us_lin, J_uu_lin];


% ****************** RHS ************************************** %

f_u_vec = zeros(Nx_pts*Nt_pts,1);
f_sigma_vec = zeros(Nx_pts*Nt_pts,1);

f = [f_sigma_vec; f_u_vec];
 
rhs_sigma = f_sigma_vec;
rhs_u = f_u_vec;

[rhs_sym, Hess_J_sym] = apply_bdy_cond(Nx_pts, Nt_pts, Hess_J, f, bdy_cond, u0, bdy_left, bdy_right, 1);

% *********** EIGENVALUES w/ & w/o BOUNDARY CONDITIONS ********* %

[V_sym, D_sym] = eig(Hess_J_sym);
D_sym = sort(diag(D_sym));

% [V, D] = eig(Hess_J);
% D = sort(diag(D));

% ****************** PLOT ************************************* %

font_size = 20; 

figure;
scatter(1:length(D_sym), D_sym, '.', 'b');
ylabel('eigenvalues', 'interpreter','latex', 'FontSize', font_size);

% title('Eigenvalues Hess_J sym, 11 x 11');

fprintf('condition number symmetrised hessian %d\n', cond(Hess_J_sym));
% condition number of the non symmetrised case so much higher, 
% because matrix probably more or less singular, many eigenvalues
% VERY close to zero
% figure;
% scatter(1:length(D), D, '.');
% % title('Eigenvalues Hess_J sym, 11 x 11');
% fprintf('condition number hessian w/o bdy cond %d\n', cond(Hess_J));



