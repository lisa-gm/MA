% ****************************************************************** %
% ******************* MAKING PRETTY GRAPHS ************************* %
% ****************************************************************** %

% SHOWING PICTURES OF INTERPOLATION & RECSTRICTION OPERATOR

clear all;
close all;

% CHOOSE GENERAL PARAMTERS: 
levels = 2;

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
Nx_elem_list(1) = 8;

% time
T = 1;
Nt_elem_list(1) = 8;


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

x_vec_c = 0:hx_list(2):S;
t_vec_c = 0:ht_list(2):T;

Nx_pts_list = Nx_elem_list + 1;
Nt_pts_list = Nt_elem_list + 1;

tot_pts = Nx_pts_list(1)*Nt_pts_list(1);
tot_pts_coarse = Nx_pts_list(2)*Nt_pts_list(2);

% ************************ BOUNDARY CONDITIONS ************************** %

%bdy_cond = 'Dirichlet';
bdy_cond = 'Neumann';

% *************** INITIALISE U & SIGMA w/ bdy * FINE ******************** %

%u0 = sin(pi*x_vec);
u0 = transpose(max(1-2*x_vec_f, 0));

bdy_left = zeros(Nt_pts_list(1),1);
bdy_right = zeros(Nt_pts_list(1),1);

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

% ************** INITIALISE U & SIGMA w/ bdy * COARSE ******************* %

%u0 = sin(pi*x_vec);
u0_coarse = transpose(max(1-2*x_vec_c, 0));

bdy_left_coarse = zeros(Nt_pts_list(2),1);
bdy_right_coarse = zeros(Nt_pts_list(2),1);

sigma_bdy_mat_coarse = 0*ones(Nx_pts_list(2), Nt_pts_list(2));
u_bdy_mat_coarse = zeros(Nx_pts_list(2), Nt_pts_list(2));

for ts=1:Nt_pts_list(2)
    u_bdy_mat_coarse(:,ts) = u0_coarse;
end

if(strcmp(bdy_cond, 'Dirichlet'))
    u_bdy_mat_coarse(1, :) = bdy_left_coarse;
    u_bdy_mat_coarse(end, :) = bdy_right_coarse;
end 

if(strcmp(bdy_cond, 'Neumann'))
    sigma_bdy_mat_coarse(1, :) = bdy_left_coarse;
    sigma_bdy_mat_coarse(end, :) = bdy_right_coarse;
end

sigma_bdy_coarse = reshape(sigma_bdy_mat_coarse, [tot_pts_coarse,1]);
u_bdy_coarse = reshape(u_bdy_mat_coarse, [tot_pts_coarse,1]);

sol_init_coarse = [sigma_bdy_coarse; u_bdy_coarse];

% ****************** CONSTRUCTING MATRICES * FINE ****************** %

[J_ss_lin, J_su_lin, J_us_lin, J_uu_lin] = lin_part_mat_no_bdy_lsfem (S, Nx_elem_list(1), T, Nt_elem_list(1), c1 ,c2, diff_const);
Hess_J = [J_ss_lin, J_su_lin; J_us_lin, J_uu_lin];

% ****************** CONSTRUCTING MATRICES * COARSE ****************** %

[J_ss_lin_coarse, J_su_lin_coarse, J_us_lin_coarse, J_uu_lin_coarse] = lin_part_mat_no_bdy_lsfem (S, Nx_elem_list(2), T, Nt_elem_list(2), c1 ,c2, diff_const);
Hess_J_coarse = [J_ss_lin_coarse, J_su_lin_coarse; J_us_lin_coarse, J_uu_lin_coarse];


% ****************** RHS * FINE ************************************** %

f_u_mat = zeros(Nx_pts_list(1), Nt_pts_list(1));
f_u_vec = reshape(f_u_mat, [Nx_pts_list(1)*Nt_pts_list(1),1]);

f_sigma_vec = zeros(Nx_pts_list(1)*Nt_pts_list(1),1);

f = [f_sigma_vec; f_u_vec];
 
rhs_sigma = f_sigma_vec;
rhs_u = f_u_vec;

[rhs_sym, Hess_J_sym] = apply_bdy_cond(Nx_pts_list(1), Nt_pts_list(1), Hess_J, f, bdy_cond, u0, bdy_left, bdy_right, 1);

% ****************** RHS * COARSE ************************************** %

f_u_vec_coarse = zeros(Nx_pts_list(2)*Nt_pts_list(2),1);
f_sigma_vec_coarse = zeros(Nx_pts_list(2)*Nt_pts_list(2),1);

f_coarse = [f_sigma_vec_coarse; f_u_vec_coarse];

[rhs_sym_coarse, Hess_J_sym_coarse] = apply_bdy_cond(Nx_pts_list(2), Nt_pts_list(2), Hess_J_coarse, f_coarse, bdy_cond, u0_coarse, bdy_left_coarse, bdy_right_coarse, 1);

% ************* DIRECT SOLVE * FINE ******************* %

sol_direct = Hess_J_sym \ rhs_sym;

sigma_direct = sol_direct(1:tot_pts);
u_direct = sol_direct(tot_pts+1:end);

% ************* DIRECT SOLVE * COARSE ******************* %

sol_direct_coarse = Hess_J_sym_coarse \ rhs_sym_coarse;

sigma_direct_coarse = sol_direct_coarse(1:tot_pts_coarse);
u_direct_coarse = sol_direct_coarse(tot_pts_coarse+1:end);

% get interpolation matrices
[I, R] = set_up_interpolation_op_SP_TIME_no_bdy(Nx_elem_list, Nt_elem_list, levels);

% ************************ RESTRICT ******************************** %
sol_rest_coarse = R{1}*sol_direct;

sigma_rest_coarse = sol_rest_coarse(1:tot_pts_coarse);
u_rest_coarse = sol_rest_coarse(tot_pts_coarse+1:end);

% ************************ INPTERPOLATE **************************** %
sol_int_coarse = I{1}*sol_direct_coarse;

sigma_int_coarse = sol_int_coarse(1:tot_pts);
u_int_coarse = sol_int_coarse(tot_pts+1:end);

% ************************ PLOTTING ******************************** %

% reshape
u_mat_direct = zeros(Nx_pts_list(1), Nt_pts_list(1));
sigma_mat_direct = zeros(Nx_pts_list(1), Nt_pts_list(1));

u_mat_int = zeros(Nx_pts_list(1), Nt_pts_list(1));
sigma_mat_int = zeros(Nx_pts_list(1), Nt_pts_list(1));

for ts=1:Nt_pts_list(1)
    sigma_mat_direct(:,ts) = sigma_direct((ts-1)*Nx_pts_list(1)+1:ts*Nx_pts_list(1));
    u_mat_direct(:, ts) = u_direct((ts-1)*Nx_pts_list(1)+1:ts*Nx_pts_list(1));
    
    sigma_mat_int(:,ts) = sigma_int_coarse((ts-1)*Nx_pts_list(1)+1:ts*Nx_pts_list(1));
    u_mat_int(:, ts) = u_int_coarse((ts-1)*Nx_pts_list(1)+1:ts*Nx_pts_list(1));    
end

% restriction or interpolation
u_mat_rest_coarse = zeros(Nx_pts_list(2), Nt_pts_list(2));
sigma_mat_rest_coarse = zeros(Nx_pts_list(2), Nt_pts_list(2));

u_mat_direct_coarse = zeros(Nx_pts_list(2), Nt_pts_list(2));
sigma_mat_direct_coarse = zeros(Nx_pts_list(2), Nt_pts_list(2));

for ts=1:Nt_pts_list(2)
    sigma_mat_direct_coarse(:,ts) = sigma_direct_coarse((ts-1)*Nx_pts_list(2)+1:ts*Nx_pts_list(2));
    u_mat_direct_coarse(:, ts) = u_direct_coarse((ts-1)*Nx_pts_list(2)+1:ts*Nx_pts_list(2));    
    
    sigma_mat_rest_coarse(:,ts) = sigma_rest_coarse((ts-1)*Nx_pts_list(2)+1:ts*Nx_pts_list(2));
    u_mat_rest_coarse(:, ts) = u_rest_coarse((ts-1)*Nx_pts_list(2)+1:ts*Nx_pts_list(2));
end

font_size = 15;

figure;
subplot(1,2,1);
mesh(x_vec_c, t_vec_c, transpose(u_mat_direct_coarse));
xlabel('space','FontSize', font_size);
ylabel('time','FontSize', font_size);
zlabel('u','FontSize', font_size);
colormap winter;
%zlim([min(u), max(u)]);

subplot(1,2,2);
mesh(x_vec_f, t_vec_f, transpose(u_mat_int)); 
xlabel('space','FontSize', font_size);
ylabel('time','FontSize', font_size);
zlabel('u','FontSize', font_size);
colormap winter;
%zlim([min(u), max(u)]);


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







