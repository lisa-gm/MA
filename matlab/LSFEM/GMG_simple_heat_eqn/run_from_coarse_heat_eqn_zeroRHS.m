% ****************************************************************** %
% *************** SOLVING LSFEM HEAT EQN WITH F = 0 **************** %
% ******************** USING GEOMETRIC MULTIGRID ******************* %
% ****************************************************************** %

% SETTING UP GENERAL PARAMETERS + MULTIGRID

%clear all;
close all;

% CHOOSE GENERAL PARAMTERS: 

% MG parameters
max_iter = 400;                          % number of newton iterations
eps = 10^(-9);                           % tolerance
levels = 2; 
%smoother = 'GaussSeidel';
smoother = 'Jacobi_LS_extended_SP_T';
%smoother = 'two_block_GaussSeidel';
%smoother = 'line_smoother';

% eqn parameters
diff_const = 1;                          % diffusion constant, define through sigma = diff_const*nabla(u) 

% LSFEM parameters
c1 = 1;                                 % constants in front of J = c1*|| ....-f(u) ||^2 + c2*|| ...||^2
c2 = 1;                                 % need to multiply with them in computation of each integral 
                                        
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
Nx_elem_list(1) = 8;

% time
T = 1;
Nt_elem_list(1) = 8;


% adaptively refine mesh by factor of 2, add in pts in space & time
for j=1:levels
    hx_list(j) = S/Nx_elem_list(j);
    ht_list(j) = T/Nt_elem_list(j);
    
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

bdy_left =  zeros(Nt_pts_list(1),1);
bdy_right = zeros(Nt_pts_list(1),1);

% ***************** INITIALISE U & SIGMA w/ bdy ******************************* %

sigma_bdy_mat = 0*ones(Nx_pts_list(1), Nt_pts_list(1));
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

% ********************* CONSTRUCTING MATRICES *********************** %

[J_ss_lin, J_su_lin, J_us_lin, J_uu_lin] = lin_part_mat_no_bdy_lsfem (S, Nx_elem_list(1), T, Nt_elem_list(1), c1 ,c2, diff_const);
Hess_J = [J_ss_lin, J_su_lin; J_us_lin, J_uu_lin];
% 

% ****************** RHS ******************************************** %

f_u_mat = zeros(Nx_pts_list(1), Nt_pts_list(1));
f_u_vec = reshape(f_u_mat, [Nx_pts_list(1)*Nt_pts_list(1),1]);

f_sigma_vec = zeros(Nx_pts_list(1)*Nt_pts_list(1),1);

f = [f_sigma_vec; f_u_vec];
 
rhs_sigma = f_sigma_vec;
rhs_u = f_u_vec;

[rhs_sym, Hess_J_sym] = apply_bdy_cond(Nx_pts_list(1), Nt_pts_list(1), Hess_J, f, bdy_cond, u0, bdy_left, bdy_right, 1);

% ********************** MULTIGRID ********************************** %

%[sol_mg, it_sol_mg, it_res_mg] = V_cycle_plain(Nx_elem_list, hx_list, Nt_elem_list, ht_list, f, sol_init, c1, c2, diff_const, bdy_cond, u0, bdy_left, bdy_right, levels, max_iter, smoother, eps);
[sol_mg, it_sol_mg, it_res_mg] = V_cycle_galerkin_assembly(Nx_elem_list, hx_list, Nt_elem_list, ht_list, Hess_J, f, sol_init, bdy_cond, u0, bdy_left, bdy_right, levels, max_iter, smoother, eps);

file_name = ['it_res_mg_',num2str(Nx_elem_list(1)),'by', num2str(Nt_elem_list(1)),'_fine_elem_c2_', num2str(c2),'.mat'];
save(file_name, 'it_res_mg');

sigma_mg = sol_mg(1:tot_pts);
u_mg = sol_mg(tot_pts+1:end);

% ************* DIRECT SOLVE FOR COMPARISON ************************* %

sol_direct = Hess_J_sym \ rhs_sym;

sigma_direct = sol_direct(1:tot_pts);
u_direct = sol_direct(tot_pts+1:end);

% **************** COMPARISON MG AND DIRECT SOLVER ***************** %

norm_diff_sol = norm(sol_direct - sol_mg);
norm_diff_u = norm(u_direct - u_mg);

fprintf('\n');
fprintf('final norm residual after %d iterations: %d\n', length(it_res_mg)-1, it_res_mg(end));
fprintf('difference MG vs DIRECT solution : %d\n', norm_diff_sol);
%fprintf('\n');

% compute convergence coefficient : || u - u_k || / || u - u_k-1 ||
conv_coeff = zeros(1,size(it_sol_mg,2)-1);
conv_coeff_res = zeros(1, length(it_res_mg)-1);

for iter=1:length(it_res_mg)-1
        conv_coeff(iter) = norm(sol_direct - it_sol_mg(:,iter+1))/norm(sol_direct - it_sol_mg(:,iter));
        conv_coeff_res(iter) = it_res_mg(iter+1)/it_res_mg(iter);    
end

file_name = ['conv_coeff_res_',num2str(Nx_elem_list(1)),'by', num2str(Nt_elem_list(1)),'_fine_elem_c2_', num2str(c2),'.mat'];
save(file_name, 'conv_coeff_res');


% figure;
% plot(1:length(conv_coeff), conv_coeff);
% xlabel('iterations')
% ylim([0.5,1.1]);
% title('convergence coefficient: || u - u_k || / || u - u_{k-1} ||');

font_size = 20;

figure;
semilogy(1:length(it_res_mg), it_res_mg, 'b', 'LineWidth', 2);
xlim([0, length(it_res_mg)]);
xlabel('iterations', 'interpreter','latex', 'FontSize', font_size);
ylabel('residual norm', 'interpreter','latex', 'FontSize', font_size);
%title('norm residual, log plot', 'interpreter','latex', 'FontSize', font_size);

figure;
plot(1:length(conv_coeff_res),conv_coeff_res, 'b', 'LineWidth', 2);
xlim([0, length(conv_coeff_res)]);
xlabel('k', 'interpreter','latex', 'FontSize', font_size);
ylabel('convergence coefficient', 'interpreter','latex', 'FontSize', font_size);
%title('convergence coefficient residual : || res_k || / || res_{k-1} ||', 'interpreter','latex', 'FontSize', font_size);

fprintf('average value res conv coeff : %d \n', mean(conv_coeff_res));
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
subplot(1,2,1);
mesh(x_vec_f, t_vec_f, transpose( u_bdy_mat));
xlabel('space', 'interpreter','latex', 'FontSize', font_size);
ylabel('time', 'interpreter','latex', 'FontSize', font_size);
zlabel('u', 'interpreter','latex', 'FontSize', font_size);
%zlim([min(u), max(u)]);
% title('u initial guess', 'interpreter','latex', 'FontSize', font_size);

subplot(1,2,2);
mesh(x_vec_f, t_vec_f, transpose(u_mat_mg));
xlabel('space', 'interpreter','latex', 'FontSize', font_size);
ylabel('time', 'interpreter','latex', 'FontSize', font_size);
zlabel('u', 'interpreter','latex', 'FontSize', font_size);
%zlim([min(u), max(u)]);
%title('heat equation LSFEM, f=0, MULTIGRID V-CYCLE', 'interpreter','latex', 'FontSize', font_size);

% subplot(1,3,3);
% mesh(x_vec_f, t_vec_f, transpose(u_mat_direct)); 
% xlabel('space', 'interpreter','latex', 'FontSize', font_size);
% ylabel('time', 'interpreter','latex', 'FontSize', font_size);
% zlabel('u', 'interpreter','latex', 'FontSize', font_size);
% %zlim([min(u), max(u)]);
% %title('u direct', 'interpreter','latex', 'FontSize', font_size);


figure;
subplot(1,3,1);
mesh(x_vec_f, t_vec_f, transpose(sigma_bdy_mat)); 
xlabel('space', 'interpreter','latex', 'FontSize', font_size);
ylabel('time', 'interpreter','latex', 'FontSize', font_size);
zlabel('sigma', 'interpreter','latex', 'FontSize', font_size);
%zlim([min(u), max(u)]);
%title('SIGMA,  sigma direct', 'interpreter','latex', 'FontSize', font_size);

subplot(1,3,2);
mesh(x_vec_f, t_vec_f, transpose(sigma_mat_mg)); 
xlabel('space', 'interpreter','latex', 'FontSize', font_size);
ylabel('time', 'interpreter','latex', 'FontSize', font_size);
zlabel('sigma', 'interpreter','latex', 'FontSize', font_size);
%zlim([min(u), max(u)]);
%title('SIGMA, sigma mg', 'interpreter','latex', 'FontSize', font_size);

subplot(1,3,3);
mesh(x_vec_f, t_vec_f, transpose(sigma_mat_direct)); 
xlabel('space', 'interpreter','latex', 'FontSize', font_size);
ylabel('time', 'interpreter','latex', 'FontSize', font_size);
zlabel('sigma', 'interpreter','latex', 'FontSize', font_size);
%zlim([min(u), max(u)]);
%title('SIGMA, direct', 'interpreter','latex', 'FontSize', font_size);


% figure;
% subplot(1,2,1);
% mesh(x_vec_f, t_vec_f, transpose(sigma_mat_direct-sigma_mat_mg)); 
% xlabel('space');
% ylabel('time');
% zlabel('diff sigma');
% %zlim([min(u), max(u)]);
% title('sigma direct - sigma mg');
% 
% subplot(1,2,2);
% mesh(x_vec_f, t_vec_f, transpose(u_mat_direct-u_mat_mg)); 
% xlabel('space');
% ylabel('time');
% zlabel('diff u');
% %zlim([min(u), max(u)]);
% title('u direct - u mg');

