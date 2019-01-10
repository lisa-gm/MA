% *********************************************************************** %
% ************** EIGENFUNCTION DECOMPOSITION ERROR ********************** %
% *********************************************************************** %

% find out why convergence rate is so low, which eigenmodes don't get
% reduced well? therefore:

% let e = u - u_k where u_k : k-th iterate
% let v_1, ...,v_n be the eigenvectors of Hess_J, then
% there exists unique sum mu_j * v_j = e
% now find coefficients mu_j, check which ones are large, correspond to
% eigennodes that don't get reduced well, can also check which eigenvalue
% they belong to

% ie. e = mu_1*v_1 + ... + mu_n*v_n
% ---> mu = V \ e where V has columns consisting of eigenvectors

% *********************************************************************** %

% SETTING UP GENERAL PARAMETERS + MULTIGRID

clear all;
close all;

% CHOOSE GENERAL PARAMTERS: 

% MG parameters
max_iter = 100;                          % number of newton iterations
eps = 10^(-9);                           % tolerance
levels = 2; 
smoother = 'GaussSeidel_LS';
%smoother = 'Jacobi_LS_extended_SP_T';
%block_size = 1;
omega = 0.6656;

% eqn parameters
diff_const = 0.1;                          % diffusion constant, define through sigma = diff_const*nabla(u) 

% LSFEM parameters
c1 = 2;                                 % constants in front of J = c1*|| ....-f(u) ||^2 + c2*|| ...||^2
c2 = 2;                                 % need to multiply with them in computation of each integral 
                                        
% ********** SETTING UP DOMAIN PARAMETERS ***************************** %
% from fine to coarse
Nx_elem_list = zeros(1,levels);
Nt_elem_list = zeros(1,levels);

hx_list = zeros(1,levels);
ht_list = zeros(1,levels);

% start with COARSEST LEVEL, then adaptively refine to have nested meshes !! 
% we already know number of total levels, but 

% space 
S = 1;
Nx_elem_list(1) = 10;

% time
T = 1;
Nt_elem_list(1) = 10;


% adaptively refine mesh by factor of 2, add in pts in space & time
for j=1:levels
    hx_list(j) = S/Nx_elem_list(j);
    ht_list(j) = S/Nt_elem_list(j);
    
    if(j < levels)
    Nx_elem_list(j+1) = 2*Nx_elem_list(j);
    Nt_elem_list(j+1) = 2*Nt_elem_list(j);
    end    
end 

%%%%%%%%%%%%%%%%%%%%%%% reverse order of lists %%%%%%%%%%%%%%%%%%%%%%%%%%%
hx_list = fliplr(hx_list);
ht_list = fliplr(ht_list);

Nx_elem_list = fliplr(Nx_elem_list);
Nt_elem_list = fliplr(Nt_elem_list);
    
x_vec_f = 0:hx_list(1):S;
t_vec_f = 0:ht_list(1):T;

Nx_pts_list = Nx_elem_list + 1;
Nt_pts_list = Nt_elem_list + 1;

tot_pts = Nx_pts_list(1)*Nt_pts_list(1);

% ******** BOUNDARY CONDITIONS ON FINEST GRID ************************** %

%u0 = sin(pi*x_vec);
u0 = transpose(max(1-2*x_vec_f, 0));

%bdy_cond = 'Dirichlet';
bdy_cond = 'Neumann';

bdy_left = zeros(Nt_pts_list(1),1);
bdy_right = zeros(Nt_pts_list(1),1);

% ***************** INITIALISE U & SIGMA w/ bdy ******************************* %

sigma_bdy_mat = 0.1*ones(Nx_pts_list(1), Nt_pts_list(1));
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

% ****************** RHS ******************************************** %

f_u_mat = zeros(Nx_pts_list(1), Nt_pts_list(1));
f_u_vec = reshape(f_u_mat, [Nx_pts_list(1)*Nt_pts_list(1),1]);

f_sigma_vec = zeros(Nx_pts_list(1)*Nt_pts_list(1),1);

f = [f_sigma_vec; f_u_vec];

% ********************** MULTIGRID ********************************** %

[sol_mg, it_sol_mg, it_res_mg] = V_cycle(Nx_elem_list, hx_list, Nt_elem_list, ht_list, f, sol_init, c1, c2, diff_const, bdy_cond, u0, bdy_left, bdy_right, levels, max_iter, smoother, eps);

sigma_mg = sol_mg(1:tot_pts);
u_mg = sol_mg(tot_pts+1:end);

% ********************** DIRECT SOLVE ******************************** %

[J_ss_lin, J_su_lin, J_us_lin, J_uu_lin] = lin_part_mat_no_bdy_lsfem (S, Nx_elem_list(1), T, Nt_elem_list(1), c1 ,c2, diff_const);
Hess_J = [J_ss_lin, J_su_lin; J_us_lin, J_uu_lin];

[rhs_sym, Hess_J_sym] = apply_bdy_cond(Nx_pts_list(1), Nt_pts_list(1), Hess_J, f, bdy_cond, u0, bdy_left, bdy_right, 1);

% solve
sol_direct = Hess_J_sym \ rhs_sym;

sigma_direct = sol_direct(1:tot_pts);
u_direct = sol_direct(tot_pts+1:end);

% ********************** ANALYSIS RESULTS *************************** %

% compute convergence coefficient : || res_u_k || / || res_u_k-1 ||
conv_coeff_res = zeros(1, length(it_res_mg)-1);

for iter=1:length(it_res_mg)-1
        conv_coeff_res(iter) = it_res_mg(iter+1)/it_res_mg(iter);    
end

% figure;
% semilogy(1:length(it_res_mg), it_res_mg);
% xlabel('iterations');
% title('norm residual, log plot');
% 
% figure;
% plot(1:length(conv_coeff_res),conv_coeff_res);
% title('convergence coefficient residual : || res_k || / || res_{k-1} ||');

fprintf('num of iter: %d, final residual: %d\n', length(it_res_mg)-1, it_res_mg(end));
fprintf('avg value residual convergence coeff : %d \n', mean(conv_coeff_res));

% ********************************************************************** %
% LOOKING AT THE ERROR + EIGENVALUES + COEFFICIENTS
% compute error
e = sol_direct - sol_mg;

e_sigma = e(1:tot_pts);
e_u = e(tot_pts+1:end);

% get eigenvectors + values
[V, D] = eig(Hess_J_sym);

% determine coefficients
mu = V \ e;

% display k largest coefficients
k = 5;

[sorted_mu, ind_mu] = sort(mu, 'descend');
%fprintf('largest coefficients :\n');
coeff = transpose(sorted_mu(1:k));

% eigenvalues of k largest coefficients
%fprintf('corresponding eigenvalues :\n');
eig_val = transpose(diag(D(ind_mu(1:k), ind_mu(1:k))));

% get coefficient mu_largest*eig_vec
larg_coeff = sorted_mu(1)*V(:, ind_mu(1));

larg_coeff_sigma = larg_coeff(1:tot_pts);
larg_coeff_u = larg_coeff(tot_pts+1:end);

% 2nd largest
sec_larg_coeff = sorted_mu(2)*V(:, ind_mu(2));

sec_larg_coeff_sigma = sec_larg_coeff(1:tot_pts);
sec_larg_coeff_u = sec_larg_coeff(tot_pts+1:end);

% to double check the math
sum = zeros(2*tot_pts,1);
for i=1:2*tot_pts
    sum = sum + sorted_mu(i)*V(:, ind_mu(i));
end

norm_diff = norm(e);

% PLOT EIGENFUNCTIONS -- plot sigma and u separately ...

% reshape
e_sigma_mat = zeros(Nx_pts_list(1), Nt_pts_list(1));
e_u_mat = zeros(Nx_pts_list(1), Nt_pts_list(1));

larg_coeff_sigma_mat = zeros(Nx_pts_list(1), Nt_pts_list(1));
larg_coeff_u_mat = zeros(Nx_pts_list(1), Nt_pts_list(1));

sec_larg_coeff_sigma_mat = zeros(Nx_pts_list(1), Nt_pts_list(1));
sec_larg_coeff_u_mat = zeros(Nx_pts_list(1), Nt_pts_list(1));


for ts=1:Nt_pts_list(1)
    e_sigma_mat(:,ts) = e_sigma((ts-1)*Nx_pts_list(1)+1:ts*Nx_pts_list(1));
    e_u_mat(:, ts) = e_u((ts-1)*Nx_pts_list(1)+1:ts*Nx_pts_list(1));
    
    larg_coeff_sigma_mat(:, ts) = larg_coeff_sigma((ts-1)*Nx_pts_list(1)+1:ts*Nx_pts_list(1));
    sec_larg_coeff_sigma_mat(:, ts) = sec_larg_coeff_sigma((ts-1)*Nx_pts_list(1)+1:ts*Nx_pts_list(1));
    
    larg_coeff_u_mat(:, ts) = larg_coeff_u((ts-1)*Nx_pts_list(1)+1:ts*Nx_pts_list(1));
    sec_larg_coeff_u_mat(:, ts) = sec_larg_coeff_u((ts-1)*Nx_pts_list(1)+1:ts*Nx_pts_list(1));

end

plot_sol_and_coarse_grid(Nx_elem_list(end), Nt_elem_list(end), x_vec_f, t_vec_f, e_u_mat, e_sigma_mat);

block_size_s = 3;
block_size_t = 3;
max_iter_sm = 100;
omega = 0.6656;

fprintf('norm error before add smoothening : %d\n', norm(e));
[sol_sm, ~] = JacobiSolve_LS_extended_SP_T(Nx_pts_list(1), Nt_pts_list(1), Hess_J_sym, rhs_sym, sol_mg, max_iter_sm, omega, block_size_s, block_size_t);

e = sol_direct - sol_sm;
fprintf('norm error after add smoothening : %d\n', norm(e));

e_sigma = e(1:tot_pts);
e_u = e(tot_pts+1:end);

e_sigma_mat = zeros(Nx_pts_list(1), Nt_pts_list(1));
e_u_mat = zeros(Nx_pts_list(1), Nt_pts_list(1));


for ts=1:Nt_pts_list(1)
    e_sigma_mat(:,ts) = e_sigma((ts-1)*Nx_pts_list(1)+1:ts*Nx_pts_list(1));
    e_u_mat(:, ts) = e_u((ts-1)*Nx_pts_list(1)+1:ts*Nx_pts_list(1));
end

plot_sol_and_coarse_grid(Nx_elem_list(end), Nt_elem_list(end), x_vec_f, t_vec_f, e_u_mat, e_sigma_mat);

% figure;
% str = ['after ', num2str(length(it_res_mg)-1), ' iterations'];
% suptitle(str);
% subplot(2,3,1)
% mesh(x_vec_f, t_vec_f, transpose(e_sigma_mat));
% xlabel('space');
% ylabel('time');
% zlabel('e sigma');
% %zlim([-0.35,0.3]);
% title('error sigma');
% 
% subplot(2,3,2)
% mesh(x_vec_f, t_vec_f, transpose(larg_coeff_sigma_mat));
% xlabel('space');
% ylabel('time');
% zlabel('mu_largest*eigvec');
% %zlim([-0.35,0.3]);
% title('SIGMA, largest coeff tim eigfct');
% 
% subplot(2,3,3)
% mesh(x_vec_f, t_vec_f, transpose(sec_larg_coeff_sigma_mat));
% xlabel('space');
% ylabel('time');
% zlabel('mu_largest*eigvec');
% title('SIGMA, largest coeff tim eigfct');
% 
% subplot(2,3,4)
% mesh(x_vec_f, t_vec_f, transpose(e_u_mat));
% xlabel('space');
% ylabel('time');
% zlabel('e u');
% %zlim([-0.01,0.015]);
% title('error u');
% 
% subplot(2,3,5)
% mesh(x_vec_f, t_vec_f, transpose(larg_coeff_u_mat));
% xlabel('space');
% ylabel('time');
% zlabel('mu_largest*eigvec');
% %zlim([-0.01,0.015]);
% title('U, largest coeff tim eigfct');
% 
% subplot(2,3,6)
% mesh(x_vec_f, t_vec_f, transpose(sec_larg_coeff_u_mat));
% xlabel('space');
% ylabel('time');
% zlabel('mu_2nd largest*eigvec');
% title('U, 2nd largest coeff tim eigfct');

