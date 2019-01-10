% ********************************************************************* %
% ****************** POISSON EQN 1D W MULTIGRID *********************** %
% ********************************************************************* %

% to compare convergence rate to least squares heat eqn
clear all;
close all;

% CHOOSE GENERAL PARAMTERS: 

% MG parameters
max_iter = 200;                          % number of newton iterations
eps = 10^(-9);                           % tolerance
levels = 2; 
smoother = 'Jacobi_LS';                  

% eqn parameters
diff_const = 1;                          % diffusion constant, define through sigma = diff_const*nabla(u) 

% LSFEM parameters
c1 = 2;                                 % constants in front of J = c1*|| ....-f(u) ||^2 + c2*|| ...||^2
c2 = 2;                                 % need to multiply with them in computation of each integral 
            
% ********** SETTING UP DOMAIN PARAMETERS ***************************** %
% from fine to coarse
Nx_elem_list = zeros(1,levels);
hx_list = zeros(1,levels);

% start with COARSEST LEVEL, then adaptively refine to have nested meshes !! 
% we already know number of total levels, but 

% space 
S = 1;
Nx_elem_list(1) = 10;

% adaptively refine mesh by factor of 2, add in pts in space & time
for j=1:levels
    hx_list(j) = S/Nx_elem_list(j);
    
    if(j < levels)
    Nx_elem_list(j+1) = 2*Nx_elem_list(j);
    end    
end 

%%%%%%%%%%%%%%%%%%%%%%% reverse order of lists %%%%%%%%%%%%%%%%%%%%%%%%%%%

hx_list = fliplr(hx_list);
Nx_elem_list = fliplr(Nx_elem_list);   
x_vec_f = 0:hx_list(1):S;

Nx_pts_list = Nx_elem_list + 1;

% ******** BOUNDARY CONDITIONS ON FINEST GRID ************************** %

bdy_left = 1;
bdy_right = 1;

% ***************** INITIALISE U & SIGMA w/ bdy ******************************* %
sigma = zeros(Nx_pts_list(1),1);
u = zeros(Nx_pts_list(1),1); 
u(1) = bdy_left;
u(2) = bdy_right;

sol_init = [sigma; u];

% ********************* CONSTRUCTING MATRICES *********************** %

% same as reaction diffusion, just took out u_t part
[J_ss_lin, J_su_lin, J_us_lin, J_uu_lin] = poisson_mat_no_bdy_lsfem (S, Nx_elem_list(1), diff_const);
Hess_J = [J_ss_lin, J_su_lin; J_us_lin, J_uu_lin];


% ****************** RHS ******************************************** %

f_sigma_vec = zeros(Nx_pts_list(1),1);
%f_u_vec = 2*ones(Nx_pts_list(1),1);
f_u_vec = transpose(x_vec_f);
f = [f_sigma_vec; f_u_vec];

rhs = poisson_rhs_no_bdy_lsfem (S, Nx_elem_list(1), f);

rhs(Nx_pts_list(1)+1) = bdy_left;
rhs(end) = bdy_right;

[rhs_sym, Hess_J_sym] = apply_bdy_cond_poisson(Nx_pts_list(1), Hess_J, rhs, bdy_left, bdy_right, 1);

% ********************** JACOBI ITER ********************************** %

[sol_Jac, sol] = JacobiSolve_LS(Hess_J_sym, rhs_sym, sol_init, 300);

sigma_Jac = sol_Jac(1:Nx_pts_list(1));
u_Jac = sol_Jac(Nx_pts_list(1)+1:end);

% ********************** MULTIGRID ********************************** %

[sol_mg, it_sol_mg, it_res_mg] = V_cycle_plain_poisson(Nx_elem_list, hx_list, f, sol_init, diff_const, bdy_left, bdy_right, levels, max_iter, smoother, eps);
sol_GS = two_block_GaussSeidel(Hess_J_sym, rhs_sym, sol_init, max_iter, eps);

sigma_GS = sol_GS(1:Nx_pts_list(1));
u_GS = sol_GS(Nx_pts_list(1)+1:end);

sigma_mg = sol_mg(1:Nx_pts_list(1));
u_mg = sol_mg(Nx_pts_list(1)+1:end);

% ************* DIRECT SOLVE FOR COMPARISON ************************* %

sol_direct = Hess_J_sym \ rhs_sym;

sigma_direct = sol_direct(1:Nx_pts_list(1));
u_direct = sol_direct(Nx_pts_list(1)+1:end);

% **************** COMPARISON MG AND DIRECT SOLVER ***************** %

% norm_diff_sol = norm(sol_direct - sol_mg);
% norm_diff_u = norm(u_direct - u_mg);
% 
% fprintf('\n');
% fprintf('final norm residual after %d iterations: %d\n', length(it_res_mg)-1, it_res_mg(end));
% fprintf('difference MG vs DIRECT solution : %d\n', norm_diff_sol);
% %fprintf('\n');


conv_coeff_res = zeros(1, length(it_res_mg)-1);

for iter=1:length(it_res_mg)-1
        conv_coeff_res(iter) = it_res_mg(iter+1)/it_res_mg(iter);    
end

% figure;
% semilogy(1:length(it_res_mg), it_res_mg);
% xlabel('iterations');
% title('norm residual, log plot');
% 
figure;
plot(1:length(conv_coeff_res),conv_coeff_res);
title('convergence coefficient residual : || res_k || / || res_{k-1} ||');

fprintf('average value res conv coeff : %d \n', mean(conv_coeff_res));
% ************************ PLOTTING ******************************** %

figure;
subplot(1,2,1);
plot(x_vec_f, sigma_direct);
title('1D Laplace LSFEM, direct solve, SIGMA');
subplot(1,2,2);
plot(x_vec_f, u_direct);
title('1D Laplace LSFEM, direct solve');


figure;
subplot(1,2,1);
plot(x_vec_f, sigma_GS);
title('1D Laplace LSFEM, 2 block Gauss Seidel, SIGMA');
subplot(1,2,2);
plot(x_vec_f, u_GS);
title('1D Laplace LSFEM, 2 block Gauss Seidel');

figure;
subplot(1,2,1);
plot(x_vec_f, sigma_mg);
title('1D Laplace LSFEM, multigrid, SIGMA');
subplot(1,2,2);
plot(x_vec_f, u_mg);
title('1D Laplace LSFEM, multigrid');




