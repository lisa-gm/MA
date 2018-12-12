% ****************************************************************** %
% *************** SOLVING LSFEM HEAT EQN WITH F = 0 **************** %
% ******************** USING GEOMETRIC MULTIGRID ******************* %
% ****************************************************************** %

% SETTING UP GENERAL PARAMETERS + MULTIGRID

clear all;
close all;

% CHOOSE GENERAL PARAMTERS: 

% MG parameters
max_iter = 100;                          % number of newton iterations
eps = 10^(-9);                           % tolerance
levels = 2; 
smoother = 'GaussSeidel';
%smoother = 'Jacobi_LS_extended';

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

% get inner indices
% [inner_ind_sigma, inner_ind_u] = inner_ind(bdy_cond, Nx_elem, Nt_elem);

% ********************* CONSTRUCTING MATRICES *********************** %

[J_ss_lin, J_su_lin, J_us_lin, J_uu_lin] = lin_part_mat_no_bdy_lsfem (S, Nx_elem_list(1), T, Nt_elem_list(1), c1 ,c2, diff_const);
Hess_J = [J_ss_lin, J_su_lin; J_us_lin, J_uu_lin];
% 

% ****************** RHS ******************************************** %

f_u_mat = zeros(Nx_pts_list(1), Nt_pts_list(1));
f_u_vec = reshape(f_u_mat, [Nx_pts_list(1)*Nt_pts_list(1),1]);

f_sigma_vec = zeros(Nx_pts_list(1)*Nt_pts_list(1),1);

f = [f_sigma_vec; f_u_vec];

% 
% % compute update for RHS
% 
% Hess_J_inner_big = zeros(2*tot_pts);
% Hess_J_inner_big([inner_ind_sigma, tot_pts+inner_ind_u], [inner_ind_sigma, tot_pts+inner_ind_u]) = Hess_J_inner;
% 
% %bdy_contrib = (Hess_J - Hess_J_inner_big)*[sigma_bdy; u_bdy];
% 
rhs_sigma = f_sigma_vec;
rhs_u = f_u_vec;

[rhs_sym, Hess_J_sym] = apply_bdy_cond(Nx_pts_list(1), Nt_pts_list(1), Hess_J, f, bdy_cond, u0, bdy_left, bdy_right, 1);

% ************************ DAMPED JACOBI **************************** %
% computing optimal damping parameter for blocked Jacobi
% omega_opt = 2 / (lambda_min(D^-1*Hess_J) + lambda_max(D^-1*Hess_J))

% P = zeros(2*tot_pts);
% % compute local inverses and write them into the right places
% for k=1:tot_pts
%     loc_inv = inv(Hess_J([k, k+tot_pts],[k, k+tot_pts]));
%     P([k, k+tot_pts],[k, k+tot_pts]) = loc_inv;
% end
% 
% % compute eigenvalues of P*Hess_J
% [~, D] = eig(P*Hess_J);
% D = sort(diag(D));
% omega = 2/( D(end)); % =0.6656

omega = 0.6656;
% ********************** MULTIGRID ********************************** %

%sigma = sigma_bdy;
%u = u_bdy;

%sol_mg = GaussSeidelSolve_LS(Hess_J_sym, rhs_sym, sol_init, max_iter);
%[sol_mg, it_sol_mg] = JacobiSolve_LS(Hess_J_sym, rhs_sym, sol_init, max_iter);

block_size = 1;
%[sol_mg, it_sol_mg] = JacobiSolve_LS_extended(Nx_pts_list(1), Nt_pts_list(1), Hess_J_sym, rhs_sym, sol_init, max_iter, omega, block_size);


%sol_mg_inner = V_cycle(Hess_J_inner, rhs_inner, [], levels, max_iter, smoother, eps);

%moother = 'Jacobi_LS';
%[sol_mg, it_sol_mg] = V_cycle(Nx_elem_list, hx_list, Nt_elem_list, ht_list, f, sol_init, c1, c2, diff_const, bdy_cond, u0, bdy_left, bdy_right, levels, max_iter, smoother, eps);

%smoother = 'GaussSeidel_LS';
[sol_mg, it_sol_mg, it_res_mg] = V_cycle(Nx_elem_list, hx_list, Nt_elem_list, ht_list, f, sol_init, c1, c2, diff_const, bdy_cond, u0, bdy_left, bdy_right, levels, max_iter, smoother, eps);


%sol_mg_inner_sigma = sol_mg_inner(1:length(inner_ind_sigma));
%sol_mg_inner_u = sol_mg_inner(length(inner_ind_sigma)+1:end);

%sigma_mg = sigma_bdy;
%u_mg = u_bdy;

%sigma_mg(inner_ind_sigma) = sol_mg_inner_sigma;
%u_mg(inner_ind_u) = sol_mg_inner_u;

sigma_mg = sol_mg(1:tot_pts);
u_mg = sol_mg(tot_pts+1:end);

% ************* DIRECT SOLVE FOR COMPARISON ************************* %

sol_direct = Hess_J_sym \ rhs_sym;
%sol_direct_inner = Hess_J_inner \ rhs_inner;

%sol_direct_inner_sigma = sol_direct_inner(1:length(inner_ind_sigma));
%sol_direct_inner_u = sol_direct_inner(length(inner_ind_sigma)+1:end);

sigma_direct = sol_direct(1:tot_pts);
u_direct = sol_direct(tot_pts+1:end);

%sigma_direct(inner_ind_sigma) = sol_direct_inner_sigma;
%u_direct(inner_ind_u) = sol_direct_inner_u;

% **************** COMPARISON MG AND DIRECT SOLVER ***************** %

%norm_diff_sol = norm(sol_direct_inner - sol_mg_inner);
norm_diff_sol = norm(sol_direct - sol_mg);

norm_diff_u = norm(u_direct - u_mg);

fprintf('\n');
fprintf('difference MG vs DIRECT solution : %d\n', norm_diff_sol);
fprintf('\n');


% compute convergence coefficient : || u - u_k || / || u - u_k-1 ||
conv_coeff = zeros(1,size(it_sol_mg,2)-1);
conv_coeff_res = zeros(1, length(it_res_mg)-1);

for iter=1:length(it_res_mg)-1
        conv_coeff(iter) = norm(sol_direct - it_sol_mg(:,iter+1))/norm(sol_direct - it_sol_mg(:,iter));
        conv_coeff_res(iter) = it_res_mg(iter+1)/it_res_mg(iter);    
end

figure;
plot(1:length(conv_coeff), conv_coeff);
xlabel('iterations')
ylim([0.5,1.1]);
title('convergence coefficient: || u - u_k || / || u - u_{k-1} ||');

figure;
semilogy(1:length(it_res_mg), it_res_mg);
xlabel('iterations');
title('norm residual, log plot');

figure;
plot(1:length(conv_coeff_res),conv_coeff_res);
title('convergence coefficient residual : || res_k || / || res_{k-1} ||');

fprintf('average value convergence coeff : %d \n', mean(conv_coeff));
% ************************ PLOTTING ******************************** %

% reshape
u_mat_direct = zeros(Nx_pts_list(1), Nt_pts_list(1));
u_mat_mg = zeros(Nx_pts_list(1), Nt_pts_list(1));

sigma_mat_direct = zeros(Nx_pts_list(1), Nt_pts_list(1));
sigma_mat_mg = zeros(Nx_pts_list(1), Nt_pts_list(1));


for ts=1:Nt_pts_list(1)
    sigma_mat_direct(:,ts) = sigma_direct((ts-1)*Nx_pts_list(1)+1:ts*Nx_pts_list(1));
    u_mat_direct(:, ts) = u_direct((ts-1)*Nx_pts_list(1)+1:ts*Nx_pts_list(1));
    
    sigma_mat_mg(:, ts) = sigma_mg((ts-1)*Nx_pts_list(1)+1:ts*Nx_pts_list(1));
    u_mat_mg(:, ts) = u_mg((ts-1)*Nx_pts_list(1)+1:ts*Nx_pts_list(1));
end

figure;
subplot(1,3,1);
mesh(x_vec_f, t_vec_f, transpose( u_bdy_mat));
xlabel('space');
ylabel('time');
zlabel('u');
%zlim([min(u), max(u)]);
title('u initial guess');

subplot(1,3,2);
mesh(x_vec_f, t_vec_f, transpose( u_mat_mg));
xlabel('space');
ylabel('time');
zlabel('u');
%zlim([min(u), max(u)]);
title('heat equation LSFEM, f=0, MULTIGRID V-CYCLE');

subplot(1,3,3);
mesh(x_vec_f, t_vec_f, transpose(u_mat_direct)); 
xlabel('space');
ylabel('time');
zlabel('u');
%zlim([min(u), max(u)]);
title('u direct');


figure;
subplot(1,3,1);
mesh(x_vec_f, t_vec_f, transpose(sigma_bdy_mat)); 
xlabel('space');
ylabel('time');
zlabel('sigma');
%zlim([min(u), max(u)]);
title('SIGMA,  sigma direct');

subplot(1,3,2);
mesh(x_vec_f, t_vec_f, transpose(sigma_mat_mg)); 
xlabel('space');
ylabel('time');
zlabel('sigma');
%zlim([min(u), max(u)]);
title('SIGMA, sigma mg');

subplot(1,3,3);
mesh(x_vec_f, t_vec_f, transpose(sigma_mat_direct)); 
xlabel('space');
ylabel('time');
zlabel('sigma');
%zlim([min(u), max(u)]);
title('SIGMA, direct');


figure;
subplot(1,2,1);
mesh(x_vec_f, t_vec_f, transpose(sigma_mat_direct-sigma_mat_mg)); 
xlabel('space');
ylabel('time');
zlabel('diff sigma');
%zlim([min(u), max(u)]);
title('sigma direct - sigma mg');

subplot(1,2,2);
mesh(x_vec_f, t_vec_f, transpose(u_mat_direct-u_mat_mg)); 
xlabel('space');
ylabel('time');
zlabel('diff u');
%zlim([min(u), max(u)]);
title('u direct - u mg');