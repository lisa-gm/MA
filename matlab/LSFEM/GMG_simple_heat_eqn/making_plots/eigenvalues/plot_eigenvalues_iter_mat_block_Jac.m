% ******************************************************************* %
% *********** LOOKING AT EIGENVALUES OF ITERATION MATRIX ************ %
% ******************************************************************* %


% SETTING UP GENERAL PARAMETERS 

clear all;
%close all;

% CHOOSE GENERAL PARAMTERS: 

% eqn parameters
diff_const = 0.1;                          % diffusion constant, define through sigma = diff_const*nabla(u) 

% LSFEM parameters
c1 = 1;                                 % constants in front of J = c1*|| ....-f(u) ||^2 + c2*|| ...||^2
c2 = 100;                                 % need to multiply with them in computation of each integral 
                                        
% ********** SETTING UP DOMAIN PARAMETERS ***************************** %

test_case = 20;

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

u_bdy_mat(:,1) = u0;

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

% get inner indices
[inner_ind_sigma, inner_ind_u] = inner_ind(bdy_cond, Nx_elem, Nt_elem);

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

% ******************************************************************** %
% PARAMETERS

block_size_s = 2;
block_size_t = 2;

omega = 0.64;

max_iter = 10;


% ********************** SMOOTHERS ********************************** %

sigma = sigma_bdy;
u = u_bdy;
sol_init = [sigma; u];

% *************** SEPARATE CONSTRUCTION * SMOOTHER ****************** %

P = construct_block_Jac_mat(Nx_pts, Nt_pts, Hess_J_sym, block_size_s, block_size_t);

% sol_Jac_total = apply_smoother(Hess_J_sym, P, rhs_sym, sol_init, omega, max_iter);

% MY ITERATION MATRIX IS 

It_MAT = eye(size(P)) - omega*P*Hess_J_sym;            % ---- check its eigenvalues 


% *********** EIGENVALUES w/ & w/o BOUNDARY CONDITIONS ********* %

[V_it_mat, D_it_mat] = eig(It_MAT);
D_it_mat = sort(diag(D_it_mat));

[V, D] = eig(Hess_J_sym);
D = sort(diag(D));

% ****************** PLOT ************************************* %

font_size = 20; 

figure;
scatter(1:length(D_it_mat), D_it_mat, '.', 'b');
ylabel('eigenvalues', 'interpreter','latex', 'FontSize', font_size);


file_name = ['eigenvalues_test_case_', num2str(test_case),'.mat'];
path = '/Users/usi/Dropbox/MA/matlab/LSFEM/GMG_simple_heat_eqn/making_plots/eigenvalues/';

save([path, file_name], 'V_it_mat', 'D_it_mat', 'V', 'D');

% figure;
% scatter(1:length(D), D, '.', 'b');
% ylabel('eigenvalues Hess_J_sym', 'interpreter','latex', 'FontSize', font_size);

% title('Eigenvalues Hess_J sym, 11 x 11');

fprintf('condition number iteration matrix %d\n', cond(D_it_mat));
