% ******************************************************************** %
% ************* SET UP FOR STABILIZED GRID REFINEMENT **************** %
% ***************** + SOLVE LINEAR LSFEM, F = 0 ********************** %
% ******************************************************************** %


% SETTING UP GENERAL PARAMETERS + MULTIGRID

clear all;
close all;

% CHOOSE GENERAL PARAMTERS: 

% MG parameters
max_iter = 100;                          % number of newton iterations
eps = 10^(-9);                           % tolerance
levels = 2; 
%smoother = 'GaussSeidel';
smoother = 'Jacobi_LS_extended_SP_T';

% eqn parameters
diff_const = 1;                          % diffusion constant, define through sigma = diff_const*nabla(u) 

% LSFEM parameters
c1 = 2;                                 % constants in front of J = c1*|| ....-f(u) ||^2 + c2*|| ...||^2
c2 = 2;  

% space 
S = 1;
Nx_elem = 3;
hx = 1/3;

% time
T = 1;

% ********************************************************************** %
% SET UP NESTED GRIDS WITH 1 = diff_const*ht / hx^2
[Nx_elem_list, Nt_elem_list] = adapt_refine_mesh(hx, Nx_elem, T, diff_const, levels);

Nx_pts_list = Nx_elem_list + 1;
Nt_pts_list = Nt_elem_list + 1;

x_vec_f = 0:1/Nx_elem_list(1):S;

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

% ********************************************************************** %
[I, R] = set_up_interpolation_op_strategic(Nx_elem_list, Nt_elem_list, bdy_cond, levels);
