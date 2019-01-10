% ****************************************************************** %
% ******************* MAKING PRETTY GRAPHS ************************* %
% ****************************************************************** %

% SHOWING PICTURES OF JACOBI SMOOTHER

clear all;
close all;

% iterations block Jacobi
levels = 1;
max_iter = 9000;

eps = 10^(-9);

% BLOCK SIZE
block_size_s = 2;
block_size_t = 2;

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

% ************************ BOUNDARY CONDITIONS ************************** %

%bdy_cond = 'Dirichlet';
bdy_cond = 'Neumann';

u0 = zeros(Nx_pts_list(1),1);

bdy_left = zeros(Nx_pts_list(1),1);
bdy_right = zeros(Nx_pts_list(1),1);


% ***************** INITIALISE U & SIGMA w/ bdy ******************************* %

% create vector with random variables between [0,1] as initial guess
sigma_init_mat = rand(Nx_pts_list(1), Nt_pts_list(1));
u_init_mat = rand(Nx_pts_list(1), Nt_pts_list(1));

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

% ****************** CONSTRUCTING MATRICES * FINE ****************** %

[J_ss_lin, J_su_lin, J_us_lin, J_uu_lin] = lin_part_mat_no_bdy_lsfem (S, Nx_elem_list(1), T, Nt_elem_list(1), c1 ,c2, diff_const);
Hess_J = [J_ss_lin, J_su_lin; J_us_lin, J_uu_lin];


% ****************** RHS * FINE ************************************** %

f_u_vec = zeros(Nx_pts_list(1)*Nt_pts_list(1),1);
f_sigma_vec = zeros(Nx_pts_list(1)*Nt_pts_list(1),1);

f = [f_sigma_vec; f_u_vec];
 
rhs_sigma = f_sigma_vec;
rhs_u = f_u_vec;

[rhs_sym, Hess_J_sym] = apply_bdy_cond(Nx_pts_list(1), Nt_pts_list(1), Hess_J, f, bdy_cond, u0, bdy_left, bdy_right, 1);

% ************************ BLOCK JACOBI ***************************** %

omega = 2/3;
[~, sol_Jac] = JacobiSolve_LS_extended_SP_T(Nx_pts_list(1), Nt_pts_list(1), Hess_J_sym, rhs_sym, sol_init, max_iter, omega, block_size_s, block_size_t);

sigma_all_Jac = sol_Jac(1:tot_pts, :);
u_all_Jac = sol_Jac(tot_pts+1:end, :);

% ************************ LINE SMOOTHER **************************** %

[~, sol_lineSm] = line_smoother(Nx_pts_list(1), Nt_pts_list(1), Hess_J_sym, rhs_sym, sol_init, max_iter, eps);

sigma_all_lineSm = sol_lineSm(1:tot_pts, :);
u_all_lineSm = sol_lineSm(tot_pts+1:end, :);

% ************************ PLOTTING ******************************** %

% reshape
%u_mat_Jac_init = zeros(Nx_pts_list(1), Nt_pts_list(1));
%sigma_mat_Jac_init = zeros(Nx_pts_list(1), Nt_pts_list(1));

u_mat_Jac_3iter = zeros(Nx_pts_list(1), Nt_pts_list(1));
sigma_mat_Jac_3iter = zeros(Nx_pts_list(1), Nt_pts_list(1));

u_mat_Jac_8iter = zeros(Nx_pts_list(1), Nt_pts_list(1));
sigma_mat_Jac_8iter = zeros(Nx_pts_list(1), Nt_pts_list(1));

u_mat_Jac_20iter = zeros(Nx_pts_list(1), Nt_pts_list(1));
sigma_mat_Jac_20iter = zeros(Nx_pts_list(1), Nt_pts_list(1));

% u_mat_Jac_2000iter = zeros(Nx_pts_list(1), Nt_pts_list(1));
% sigma_mat_Jac_2000iter = zeros(Nx_pts_list(1), Nt_pts_list(1));

u_mat_Jac_init = zeros(Nx_pts_list(1), Nt_pts_list(1));
sigma_mat_Jac_init = zeros(Nx_pts_list(1), Nt_pts_list(1));

u_mat_lineSm_3iter = zeros(Nx_pts_list(1), Nt_pts_list(1));
sigma_mat_lineSm_3iter = zeros(Nx_pts_list(1), Nt_pts_list(1));

u_mat_lineSm_8iter = zeros(Nx_pts_list(1), Nt_pts_list(1));
sigma_mat_lineSm_8iter = zeros(Nx_pts_list(1), Nt_pts_list(1));

u_mat_lineSm_20iter = zeros(Nx_pts_list(1), Nt_pts_list(1));
sigma_mat_lineSm_20iter = zeros(Nx_pts_list(1), Nt_pts_list(1));

% u_mat_lineSm_2000iter = zeros(Nx_pts_list(1), Nt_pts_list(1));
% sigma_mat_lineSm_2000iter = zeros(Nx_pts_list(1), Nt_pts_list(1));

for ts=1:Nt_pts_list(1)
    sigma_mat_Jac_init(:,ts) = sigma_init((ts-1)*Nx_pts_list(1)+1:ts*Nx_pts_list(1));
    u_mat_Jac_init(:, ts) = u_init((ts-1)*Nx_pts_list(1)+1:ts*Nx_pts_list(1)); 
    
    sigma_mat_Jac_3iter(:,ts) = sigma_all_Jac((ts-1)*Nx_pts_list(1)+1:ts*Nx_pts_list(1), 3);
    u_mat_Jac_3iter(:,ts) = u_all_Jac((ts-1)*Nx_pts_list(1)+1:ts*Nx_pts_list(1), 3);
    
    sigma_mat_Jac_8iter(:,ts) = sigma_all_Jac((ts-1)*Nx_pts_list(1)+1:ts*Nx_pts_list(1), 8);
    u_mat_Jac_8iter(:,ts) = u_all_Jac((ts-1)*Nx_pts_list(1)+1:ts*Nx_pts_list(1), 8);
    
    sigma_mat_Jac_20iter(:,ts) = sigma_all_Jac((ts-1)*Nx_pts_list(1)+1:ts*Nx_pts_list(1), 20);
    u_mat_Jac_20iter(:,ts) = u_all_Jac((ts-1)*Nx_pts_list(1)+1:ts*Nx_pts_list(1), 20);
    
    %sigma_mat_Jac_2000iter(:,ts) = sigma_all((ts-1)*Nx_pts_list(1)+1:ts*Nx_pts_list(1), 2000);
    %u_mat_Jac_2000iter(:,ts) = u_all((ts-1)*Nx_pts_list(1)+1:ts*Nx_pts_list(1), 2000);
    
    sigma_mat_lineSm_3iter(:,ts) = sigma_all_lineSm((ts-1)*Nx_pts_list(1)+1:ts*Nx_pts_list(1), 3);
    u_mat_lineSm_3iter(:,ts) = u_all_lineSm((ts-1)*Nx_pts_list(1)+1:ts*Nx_pts_list(1), 3);
    
    sigma_mat_lineSm_8iter(:,ts) = sigma_all_lineSm((ts-1)*Nx_pts_list(1)+1:ts*Nx_pts_list(1), 8);
    u_mat_lineSm_8iter(:,ts) = u_all_lineSm((ts-1)*Nx_pts_list(1)+1:ts*Nx_pts_list(1), 8);
    
    sigma_mat_lineSm_20iter(:,ts) = sigma_all_lineSm((ts-1)*Nx_pts_list(1)+1:ts*Nx_pts_list(1), 20);
    u_mat_lineSm_20iter(:,ts) = u_all_lineSm((ts-1)*Nx_pts_list(1)+1:ts*Nx_pts_list(1), 20);
    
end

font_size = 15;

% figure;
% subplot(1,2,1);
% mesh(x_vec_f, t_vec_f, transpose(u_mat_Jac_init));
% xlabel('space','FontSize', font_size);
% ylabel('time','FontSize', font_size);
% zlabel('u','FontSize', font_size);
% colormap winter;
% %zlim([min(u), max(u)]);
% 
% subplot(1,2,2);
% mesh(x_vec_f, t_vec_f, transpose(u_mat_Jac_3iter)); 
% xlabel('space','FontSize', font_size);
% ylabel('time','FontSize', font_size);
% zlabel('u','FontSize', font_size);
% colormap winter;
% %zlim([min(u), max(u)]);


% figure;
% subplot(1,2,1);
% mesh(x_vec_f, t_vec_f, transpose(u_mat_Jac_init));
% xlabel('space','FontSize', font_size);
% ylabel('time','FontSize', font_size);
% zlabel('u','FontSize', font_size);
% colormap winter;
% %zlim([min(u), max(u)]);
% 
% subplot(1,2,2);
% mesh(x_vec_f, t_vec_f, transpose(u_mat_Jac_3iter)); 
% xlabel('space','FontSize', font_size);
% ylabel('time','FontSize', font_size);
% zlabel('u','FontSize', font_size);
% colormap winter;
% %zlim([min(u), max(u)]);


figure;
subplot(2,2,1);
contourf(x_vec_f, t_vec_f, transpose(u_mat_Jac_init));
caxis([-0.2,1.2]);
xlabel('space','FontSize', font_size);
ylabel('time','FontSize', font_size);
title('random initial guess','FontSize', font_size)

subplot(2,2,2);
contourf(x_vec_f, t_vec_f, transpose(u_mat_Jac_3iter)); 
caxis([-0.2,1.2]);
xlabel('space','FontSize', font_size);
ylabel('time','FontSize', font_size);
title('after 3 iterations','FontSize', font_size)

subplot(2,2,3);
contourf(x_vec_f, t_vec_f, transpose(u_mat_Jac_8iter));
xlabel('space','FontSize', font_size);
ylabel('time','FontSize', font_size);
caxis([-0.2,1.2]);
title('after 8 iterations','FontSize', font_size);

subplot(2,2,4);
contourf(x_vec_f, t_vec_f, transpose(u_mat_Jac_20iter)); 
xlabel('space','FontSize', font_size);
ylabel('time','FontSize', font_size);
caxis([-0.2,1.2]);
title('after 20 iterations','FontSize', font_size)

hp4 = get(subplot(2,2,4),'Position');
colorbar('Position', [hp4(1)+hp4(3)+0.03  hp4(2)  0.02  hp4(2)+hp4(3)*2.1]);

fprintf('norm residual Jac sm after %d iter : %d\n', length(sol_Jac(1, :)), norm(Hess_J_sym*sol_Jac(:, end) - rhs_sym));

% figure;
% subplot(1,2,1);
% mesh(x_vec_c, t_vec_c, transpose(sigma_mat_direct_coarse)); 
% xlabel('space');
% ylabel('time');
% zlabel('sigma');
% 
% subplot(1,2,2);
% mesh(x_vec_f, t_vec_f, transpose(sigma_mat_int)); 
% xlabel('space');
% ylabel('time');
% zlabel('sigma');
% %zlim([min(u), max(u)]);

figure;
subplot(2,2,1);
contourf(x_vec_f, t_vec_f, transpose(u_mat_Jac_init));
caxis([-0.2,1.2]);
xlabel('space','FontSize', font_size);
ylabel('time','FontSize', font_size);
title('random initial guess','FontSize', font_size)

subplot(2,2,2);
contourf(x_vec_f, t_vec_f, transpose(u_mat_lineSm_3iter)); 
caxis([-0.2,1.2]);
xlabel('space','FontSize', font_size);
ylabel('time','FontSize', font_size);
title('after 3 iterations','FontSize', font_size)

subplot(2,2,3);
contourf(x_vec_f, t_vec_f, transpose(u_mat_lineSm_8iter));
xlabel('space','FontSize', font_size);
ylabel('time','FontSize', font_size);
caxis([-0.2,1.2]);
title('after 8 iterations','FontSize', font_size);

subplot(2,2,4);
contourf(x_vec_f, t_vec_f, transpose(u_mat_lineSm_20iter)); 
xlabel('space','FontSize', font_size);
ylabel('time','FontSize', font_size);
caxis([-0.2,1.2]);
title('after 20 iterations','FontSize', font_size)

hp4 = get(subplot(2,2,4),'Position');
colorbar('Position', [hp4(1)+hp4(3)+0.03  hp4(2)  0.02  hp4(2)+hp4(3)*2.1]);

fprintf('norm residual line sm after %d iter : %d\n', length(sol_lineSm(1, :)), norm(Hess_J_sym*sol_lineSm(:, end) - rhs_sym));


