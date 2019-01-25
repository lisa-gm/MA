%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% LSFEM NONLINEAR RHS %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%% PARAMETERS REMOVED %%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% ADDITIONAL FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [sigma_mat, u_mat] = lsfem_nonlinear_fct(max_iter, eps, diff_const, S, Nx_elem, T, Nt_elem, ...
    bdy_cond, bdy_left, bdy_right, u0, c1, c2)

int_space = [0,S];
Nx_pts = Nx_elem+1;
hx = (int_space(2) - int_space(1))/Nx_elem;
x_vec = linspace(int_space(1), int_space(2), Nx_pts);

int_time = [0,T];
Nt_pts = Nt_elem+1;
ht = (int_time(2) - int_time(1))/Nx_elem;
t_vec = linspace(int_time(1), int_time(2), Nt_pts);

tot_pts = Nx_pts*Nt_pts;
%%%%%%%%% CONSRUCTING THE LINEAR OPERATORS %%%%%%%%%%%%%%%%%

[J_ss_lin, J_su_lin, J_us_lin, J_uu_lin] = lin_part_mat_lsfem (S, Nx_elem, T, Nt_elem, c1 , c2, diff_const);
Hess_J_linear = [J_ss_lin, J_su_lin; J_us_lin, J_uu_lin];

% *************** GETTING AN INITIAL GUESS FOR U ************************* %
sigma = 0*ones(tot_pts, 1);
u_mat = 0*ones(Nx_pts, Nt_pts);

for ts=1:Nt_pts
    u_mat(:,ts) = u0;
end
u = reshape(u_mat, [tot_pts,1]);

%[sigma, u] = get_init_guess_lsfem(S, Nx_elem, T, Nt_elem, u0, bdy_cond, bdy_left, bdy_right, a0, a1, J_ss_lin, J_su_lin, J_us_lin, J_uu_lin, c1, c2, diff_const);

if(strcmp(bdy_cond, 'Dirichlet'))
    inner_ind_u = [];
    for ind=2:Nt_pts
        inner_ind_u = [inner_ind_u, (ind-1)*Nx_pts+2:ind*Nx_pts-1];
    end
    inner_ind_sigma = 1:tot_pts;
end

if(strcmp(bdy_cond, 'Neumann'))
    inner_ind_sigma = [];
    for ind=1:Nt_pts
        inner_ind_sigma = [inner_ind_sigma, (ind-1)*Nx_pts+2:ind*Nx_pts-1];
    end
    inner_ind_u = Nx_pts+1:tot_pts;
end

%%%%%%%%%%%%%%%%%%%%%%%%% NEWTON ITERATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% damp_fact = 0;
% for iter=1:max_iter 
%    damp_fact = (1+2*damp_fact)/3;
% 
%    [sigma, u, norm_grad_J, damp_fact] = newton_step_lsfem(S, Nx_elem, T, Nt_elem, u0, bdy_cond, bdy_left, bdy_right, sigma, u, J_ss_lin, J_su_lin, J_us_lin, J_uu_lin, c1, c2, diff_const, eps, damp_fact, iter);
% 
%     if(norm_grad_J < eps)
%         fprintf('convergence criterion reached after %d iterations\n', iter);
%         break;
%     end
% end

%%%%%%%%%%%%%%%% TRUSTREGION - DOGLEG METHOD %%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up is such that we can do the multigrid call from here
max_iter_TR = 200;
TR_rad_max = 2;
TR_rad = 0.75;
eta = 0.2;

%%%%%%%%%%%%%%% MULTIGRID PARAMETERS %%%%%%%%%%%%%%%%%%
levels = 3;
smoother = 'GaussSeidel';
max_iter_mg = 20;

sigma_direct = sigma;
u_direct = u;

for iter=1:max_iter

% computing gradient and hessian    
% grad_J = grad_J_eval(hx, Nx_elem, ht, Nt_elem, bdy_cond, sigma, u, J_ss_lin, J_su_lin, J_us_lin, J_uu_lin, c1);
% Hess_J = Hess_J_eval(hx, Nx_elem, ht, Nt_elem, bdy_cond, sigma, u, J_ss_lin, J_su_lin, J_us_lin, J_uu_lin, c1);

grad_J = grad_J_eval(hx, Nx_elem, ht, Nt_elem, bdy_cond, sigma_direct, u_direct, J_ss_lin, J_su_lin, J_us_lin, J_uu_lin, c1);
Hess_J = Hess_J_eval(hx, Nx_elem, ht, Nt_elem, bdy_cond, sigma_direct, u_direct, J_ss_lin, J_su_lin, J_us_lin, J_uu_lin, c1);


grad_J_inner = grad_J([inner_ind_sigma, tot_pts+inner_ind_u]); 
Hess_J_inner = Hess_J([inner_ind_sigma, tot_pts+inner_ind_u], [inner_ind_sigma, tot_pts+inner_ind_u]);

Hess_J_linear_inner = Hess_J_linear([inner_ind_sigma, tot_pts+inner_ind_u], [inner_ind_sigma, tot_pts+inner_ind_u]);
%%%%%%%%%%%%%%%%%%%%%%%% DIRECT SOLVER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[sigma_direct, u_direct] = direct_solve_lsfem(Nx_pts, Nt_pts, inner_ind_sigma, inner_ind_u, bdy_cond, bdy_left, bdy_right, u0, Hess_J, grad_J, sigma_direct, u_direct); 
fprintf('iter: %d, norm grad J: %d\n', iter, norm(grad_J));

%%%%%%%%%%%%%%%%%%%%%%%% MULTIGRID %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check if positive definite
%[~,p] = chol(Hess_J_inner);

update = zeros(2*tot_pts,1);
update_inner = update([inner_ind_sigma, tot_pts+inner_ind_u]);


% only smoothing
% smoothing(H, -J, e, 5, smoother)
    
% geometric multigrid
%update_inner = V_cycle(Hess_J_inner, -grad_J_inner, update_inner, levels, max_iter_mg, smoother);
%update([inner_ind_sigma, tot_pts+inner_ind_u])=update_inner;

sigma = sigma + update(1:tot_pts);
u = u + update(tot_pts+1:end);

% TRUST REGION STEP
%[sigma, u, norm_grad_J, ~, TR_rad] = TR_step(TR_rad, TR_rad_max, eta, max_iter_TR, levels, smoother, max_iter_mg, S, Nx_elem, T, Nt_elem, u0, bdy_cond, bdy_left, bdy_right, sigma, u, update, grad_J, Hess_J, J_ss_lin, J_su_lin, J_us_lin, J_uu_lin, c1, c2, diff_const, iter);

    if(norm(grad_J) < eps)
        fprintf('convergence criterion reached after %d iterations\n', iter);
        break;
    end
 
end
% ************************* PLOTTING ****************************** % 

% reshape
u_mat = zeros(Nx_pts, Nt_pts);
sigma_mat = zeros(Nx_pts, Nt_pts);

u_mat_direct = zeros(Nx_pts, Nt_pts);
sigma_mat_direct = zeros(Nx_pts, Nt_pts);

for ts=1:Nt_pts
    u_mat(:, ts) = u((ts-1)*Nx_pts+1:ts*Nx_pts);
    sigma_mat(:, ts) = sigma((ts-1)*Nx_pts+1:ts*Nx_pts);
    
    u_mat_direct(:, ts) = u_direct((ts-1)*Nx_pts+1:ts*Nx_pts);
    sigma_mat_direct(:, ts) = sigma_direct((ts-1)*Nx_pts+1:ts*Nx_pts); 

end

figure
mesh(x_vec, t_vec(1:end), transpose(u_mat(:,1:end)));
xlabel('space');
ylabel('time');
zlabel('u');
%zlim([min(u), max(u)]);
title('heat equation LSFEM nonlinear MULTIGRID');
% 
% figure
% mesh(x_vec, t_vec(1:end), transpose(sigma_mat(:,1:end)));
% xlabel('space');
% ylabel('time');
% zlabel('sigma');
% %zlim([min(u), max(u)]);
% title('SIGMA heat equation LSFEM nonlinear');

figure
mesh(x_vec, t_vec(1:end), transpose(u_mat_direct(:,1:end)));
xlabel('space');
ylabel('time');
zlabel('u');
%zlim([min(u), max(u)]);
title('heat equation LSFEM nonlinear DIRECT');

end
