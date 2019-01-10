% ****************************************************************** %
% ************************ CALL LINE SMOOTHER ********************** %
% ****************************************************************** %

% ****************************************************************** %
% *************** SOLVING LSFEM HEAT EQN WITH F = 0  *************** %
% ****************************************************************** %
% SETTING UP GENERAL PARAMETERS + MULTIGRID

clear all;
%close all;

% CHOOSE GENERAL PARAMTERS: 

% MG parameters
max_iter = 4000;                            % number of iterations
eps = 10^(-8);                           % tolerance

% eqn parameters
diff_const = 0.1;                        % diffusion constant, define through sigma = diff_const*nabla(u) 

% LSFEM parameters
c1 = 2;                                 % constants in front of J = c1*|| ....-f(u) ||^2 + c2*|| ...||^2
c2 = 2;                                 % need to multiply with them in computation of each integral 
                                        
% ********** SETTING UP DOMAIN PARAMETERS ***************************** %

% space 
S = 1;
Nx_elem = 20;

Nx_pts = Nx_elem+1;
int_space = [0,S];
hx = (int_space(2) - int_space(1))/Nx_elem;
x_vec = linspace(int_space(1), int_space(2), Nx_pts);

% time
T = 1;
Nt_elem = 20;

Nt_pts = Nt_elem+1;
int_time = [0, T];
ht = (int_time(2) - int_time(1))/Nt_elem;
t_vec = linspace(int_time(1), int_time(2), Nt_pts);

tot_elem = Nx_elem*Nt_elem;
tot_pts = Nx_pts*Nt_pts;

hxht = ht*hx;

% ******************** BOUNDARY CONDITIONS ************************** %

%u0 = sin(pi*x_vec);
u0 = max(1-2*x_vec, 0);

%bdy_cond = 'Dirichlet';
bdy_cond = 'Neumann';

bdy_left = zeros(Nt_pts,1);
bdy_right = zeros(Nt_pts,1);

% ***************** INITIALISE U & SIGMA w/ bdy ******************************* %

sigma_bdy_mat = zeros(Nx_pts, Nt_pts);
u_bdy_mat = zeros(Nx_pts, Nt_pts);

for i=1:Nt_pts
    u_bdy_mat(:,1) = u0;
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

% ********************* CONSTRUCTING MATRICES *********************** %

[J_ss_lin, J_su_lin, J_us_lin, J_uu_lin] = lin_part_mat_no_bdy_lsfem (S, Nx_elem, T, Nt_elem, c1 , c2, diff_const);
Hess_J = [J_ss_lin, J_su_lin; J_us_lin, J_uu_lin];

% ****************** RHS ******************************************** %

f_u_mat = zeros(Nx_pts, Nt_pts);
f_u_vec = reshape(f_u_mat, [tot_pts,1]);

f_sigma_vec = zeros(tot_pts,1);

% compute update for RHS

rhs_sigma = f_sigma_vec;
rhs_u = f_u_vec;

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
 
% ************ MAKING THINGS SYMMETRIC ****************************** %
% get zero entries in columns where the rows are set to 1

% for each column that corresponds to a bdy entry
if(strcmp(bdy_cond, 'Dirichlet'))
    bdy_ind_sigma = [];
    bdy_ind_u = 1:Nx_pts;
    for ts=2:Nt_pts
        bdy_ind_u = [bdy_ind_u,(ts-1)*Nx_pts+1, ts*Nx_pts];
    end
end

if(strcmp(bdy_cond, 'Neumann'))
    bdy_ind_u = 1:Nx_pts;
    bdy_ind_sigma = [];
    for ts=1:Nt_pts
        bdy_ind_sigma = [bdy_ind_sigma, (ts-1)*Nx_pts+1, ts*Nx_pts];
    end
end

bdy_ind = [bdy_ind_sigma, tot_pts + bdy_ind_u];

Hess_J_sym = Hess_J;
Hess_J_sym(bdy_ind,:) = 0;

%% symmetrize the matrix

rhs_sym(:) = rhs(:) - Hess_J_sym(:,bdy_ind)*rhs(bdy_ind);
rhs_sym = rhs_sym';
Hess_J_sym(:,bdy_ind) = 0;

%% making diagonal elements 1

tempA = speye(length(Hess_J));
tempA(bdy_ind,:) = 0;
tempB = speye(length(Hess_J));
tempD = tempB-tempA;

Hess_J_sym = Hess_J_sym + tempD;

%rhs = [rhs_sigma; rhs_u] + rhs_sym;

% ********************** LINE SMOOTHER ******************************* %

sigma = sigma_bdy;
u = u_bdy;
sol_init = [sigma; u];

sol_ls = line_smoother(Nx_pts, Nt_pts, Hess_J_sym, rhs_sym, sol_init, max_iter, eps);

sigma_ls = sol_ls(1:tot_pts);
u_ls = sol_ls(tot_pts+1:end);

% ************* DIRECT SOLVE FOR COMPARISON ************************* %

sol_direct = Hess_J_sym \ rhs_sym;

sigma_direct = sol_direct(1:tot_pts);
u_direct = sol_direct(tot_pts+1:end);

% **************** COMPARISON MG AND DIRECT SOLVER ***************** %

%norm_diff_sol = norm(sol_direct_inner - sol_mg_inner);
norm_diff_sol = norm(sol_direct - [sigma_ls; u_ls]);

fprintf('\n');
fprintf('difference MG vs DIRECT solution : %d\n', norm_diff_sol);
% ************************ PLOTTING ******************************** %

% reshape
u_mat_direct = zeros(Nx_pts, Nt_pts);
u_mat_ls = zeros(Nx_pts, Nt_pts);


for ts=1:Nt_pts
    u_mat_direct(:, ts) = u_direct((ts-1)*Nx_pts+1:ts*Nx_pts);
    u_mat_ls(:, ts) = u_ls((ts-1)*Nx_pts+1:ts*Nx_pts);
end

figure
mesh(linspace(0,1,Nx_pts), linspace(0,1,Nt_pts), transpose(u_mat_ls));
xlabel('space');
ylabel('time');
zlabel('u');
%zlim([min(u), max(u)]);
title('heat equation LSFEM, f=0, LINE SMOOTHER');

figure
mesh(x_vec, t_vec, transpose(u_mat_direct));
xlabel('space');
ylabel('time');
zlabel('u');
%zlim([min(u), max(u)]);
title('heat equation LSFEM, f=0, DIRECT');

