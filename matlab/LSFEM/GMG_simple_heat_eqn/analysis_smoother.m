% ********************************************************************** %
% ************* ANALYSIS OF THE SMOOTHING PROPERTIES ******************* %
% ********************************************************************** %

% what do we have to remove in corse grid correction?

% SETTING UP GENERAL PARAMETERS + MULTIGRID

%clear all;
close all;

% CHOOSE GENERAL PARAMTERS: 

% MG parameters
max_iter = 100;                          % number of newton iterations
eps = 10^(-9);                           % tolerance
levels = 3; 
%smoother = 'ConjugateGradient';
smoother = 'Jacobi_LS';

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
% 

% space 
S = 1;
Nx_elem_list(1) = 2;

% time
T = 1;
Nt_elem_list(1) = 2;


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
%u0 = transpose(max(1-2*x_vec_f, 0));
u0 = zeros(length(Nx_pts_list(1)),1);


%bdy_cond = 'Dirichlet';
bdy_cond = 'Neumann';

bdy_left = zeros(Nt_pts_list(1),1);
bdy_right = zeros(Nt_pts_list(1),1);

% ***************** INITIALISE U & SIGMA w/ bdy ******************************* %

%sigma_bdy_mat = zeros(Nx_pts_list(1), Nt_pts_list(1));
%u_bdy_mat = zeros(Nx_pts_list(1), Nt_pts_list(1));

sigma_bdy_mat = rand(Nx_pts_list(1), Nt_pts_list(1));
u_bdy_mat = rand(Nx_pts_list(1), Nt_pts_list(1));

% only want u(x,0) = 0, rest random numbers
for ts=1:1 %Nt_pts_list(1)
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
% Hess_J_inner = Hess_J([inner_ind_sigma, tot_pts+inner_ind_u], [inner_ind_sigma, tot_pts+inner_ind_u]);

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

P = zeros(2*tot_pts);
% compute local inverses and write them into the right places
for k=1:tot_pts
    loc_inv = inv(Hess_J([k, k+tot_pts],[k, k+tot_pts]));
    P([k, k+tot_pts],[k, k+tot_pts]) = loc_inv;
end

% compute eigenvalues of P*Hess_J
[~, D] = eig(P*Hess_J);
D = sort(diag(D));
% lambda_min pretty much zero 
omega = 2/(D(1) + D(end)); % =0.6656

omega = 0.6656;
% ********************** SMOOTHER ********************************** %
block_size = 1;
Nx_pts = Nx_elem_list(1)+1;
Nt_pts = Nt_elem_list(1)+1;

%sigma = sigma_bdy;
%u = u_bdy;

%sol_GS = GaussSeidelSolve_LS(Hess_J_sym, rhs_sym, sol_init, max_iter);
[sol_J, sol_list_J] = JacobiSolve_LS(Hess_J_sym, rhs_sym, sol_init, max_iter);
%[sol_J_ex, sol_list_J_ex] = JacobiSolve_LS_extended(Nx_pts, Nt_pts, Hess_J_sym, rhs_sym, sol_init, 1, omega, block_size);


%sigma_GS = sol_GS(1:tot_pts);
%u_GS = sol_GS(tot_pts+1:end);

sigma_J = sol_J(1:tot_pts);
u_J = sol_J(tot_pts+1:end);

sol_list_u_J = sol_list_J(tot_pts+1:end, :);
sol_list_sigma_J = sol_list_J(1:tot_pts, :);


% ************ PLOTTING SOLUTION BEFORE AND AFTER ***************** %

u_mat_J_iter = [];
sigma_mat_J_iter = [];

% reshape
sigma_mat_GS = zeros(Nx_pts_list(1), Nt_pts_list(1));
sigma_mat_J = zeros(Nx_pts_list(1), Nt_pts_list(1));
sol_init_sigma_mat = zeros(Nx_pts_list(1), Nt_pts_list(1));

u_mat_GS = zeros(Nx_pts_list(1), Nt_pts_list(1));
u_mat_J = zeros(Nx_pts_list(1), Nt_pts_list(1));
sol_init_u_mat = zeros(Nx_pts_list(1), Nt_pts_list(1));

max_it = min(size(sol_list_sigma_J,2), max_iter);

for iter=1:max_it
    for ts=1:Nt_pts_list(1)
        %sigma_mat_GS(:, ts) = u_GS((ts-1)*Nx_pts_list(1)+1:ts*Nx_pts_list(1));
        sigma_mat_J(:, ts) = sol_list_sigma_J((ts-1)*Nx_pts_list(1)+1:ts*Nx_pts_list(1), iter);
        sol_init_sigma_mat(:, ts) = sigma_bdy((ts-1)*Nx_pts_list(1)+1:ts*Nx_pts_list(1));

        %u_mat_GS(:, ts) = sigma_GS((ts-1)*Nx_pts_list(1)+1:ts*Nx_pts_list(1));
        u_mat_J(:, ts) = sol_list_u_J((ts-1)*Nx_pts_list(1)+1:ts*Nx_pts_list(1), iter);
        sol_init_u_mat(:, ts) = u_bdy((ts-1)*Nx_pts_list(1)+1:ts*Nx_pts_list(1));
    end
    u_mat_J_iter{iter} = u_mat_J;
    sigma_mat_J_iter{iter} = sigma_mat_J;
    
end

% 
% % sigma
figure;
subplot(2,2,1);
mesh(x_vec_f, t_vec_f, transpose(sol_init_sigma_mat));
xlabel('space');
ylabel('time');
zlabel('sigma');
zlim([-0.2, 1.2]);
titl = 'SIGMA, heat equation LSFEM, f=0, initial sigma';
title(titl);

iter = 10;
subplot(2,2,2);
mesh(x_vec_f, t_vec_f, transpose(sigma_mat_J_iter{iter}));
xlabel('space');
ylabel('time');
zlabel('sigma');
zlim([-0.2, 1.2]);
titl = ['SIGMA heat equation LSFEM, f=0, after ', num2str(iter), ' iter of LS Jacobi sm'];
title(titl);

iter = 30;
subplot(2,2,3);
mesh(x_vec_f, t_vec_f, transpose(sigma_mat_J_iter{iter}));
xlabel('space');
ylabel('time');
zlabel('sigma');
zlim([-0.2, 1.2]);
titl = ['SIGMA heat equation LSFEM, f=0, after ', num2str(iter), ' iter of LS Jacobi sm'];
title(titl);

iter = 100;
subplot(2,2,4);
mesh(x_vec_f, t_vec_f, transpose(sigma_mat_J_iter{iter}));
xlabel('space');
ylabel('time');
zlabel('sigma');
zlim([-0.2, 1.2]);
titl = ['sigma heat equation LSFEM, f=0, after ', num2str(iter), 'iter of LS Jacobi sm'];
title(titl);


% u
figure;
subplot(2,2,1);
mesh(x_vec_f, t_vec_f, transpose(sol_init_u_mat));
xlabel('space');
ylabel('time');
zlabel('u');
zlim([-0.2, 1.2]);
titl = 'sol heat equation LSFEM, f=0, initial sigma';
title(titl);

iter = 10;
subplot(2,2,2);
mesh(x_vec_f, t_vec_f, transpose(u_mat_J_iter{iter}));
xlabel('space');
ylabel('time');
zlabel('u');
zlim([-0.2, 1.2]);
titl = ['sol heat equation LSFEM, f=0, after ', num2str(iter), ' iter of LS Jacobi sm'];
title(titl);

iter = 30;
subplot(2,2,3);
mesh(x_vec_f, t_vec_f, transpose(u_mat_J_iter{iter}));
xlabel('space');
ylabel('time');
zlabel('u');
zlim([-0.2, 1.2]);
titl = ['sol heat equation LSFEM, f=0, after ', num2str(iter), ' iter of LS Jacobi sm'];
title(titl);

iter = 100;
subplot(2,2,4);
mesh(x_vec_f, t_vec_f, transpose(u_mat_J_iter{iter}));
xlabel('space');
ylabel('time');
zlabel('u');
zlim([-0.2, 1.2]);
titl = ['sol heat equation LSFEM, f=0, after ', num2str(iter), 'iter of LS Jacobi sm'];
title(titl);

% *************** LOOKING AT THE EIGENFUNCTIONS OF ********************* %
% ********************** (I - P*HESS_J) ******************************** %

it_mat = eye(size(Hess_J)) - omega*P*Hess_J;

[V, D] = eig(it_mat); [~, perm] = sort(abs(diag(D)), 'descend');
D = D(perm, perm); V = V(:, perm);

% for now look at eigenfunctions corresponding to num largest eigenvalues\
num = 25;
max_eig_vec = V(:, 1:num);
max_eig_val = D(1:num, 1:num);

size_h = ceil(sqrt(num));
size_v = ceil(num/size_h);

% figure;
%     suptitle(['Eigenfunctions of the ', num2str(num), ' largest eigenvalues']);
%     for i=1:size_h
%         for j=1:size_v
%             if((i-1)*size_v+j <= num)
%                 
%             subplot(size_v, size_h, (i-1)*size_v+j)
%             V_sigma = reshape(max_eig_vec(1:tot_pts, (i-1)*size_v+j), [Nx_pts_list(1), Nt_pts_list(1)]);
%             V_u = reshape(max_eig_vec(tot_pts+1:end, (i-1)*size_v+j), [Nx_pts_list(1), Nt_pts_list(1)]);
%             mesh(1:2*Nx_pts_list(1), 1:Nt_pts_list(1), [V_sigma', V_u']);
%             %zlabel('V(:,1)');
%             zlim([-0.2,0.2]);
%             str = ['EigVal ', num2str(max_eig_val((i-1)*size_v+j, (i-1)*size_v+j))];
%             title(str);
%             end
%         end
%     end
