%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% LSFEM NONLINEAR RHS %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%% PARAMETERS REMOVED %%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% ADDITIONAL FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [sigma_mat, u_mat, conv_data] = lsfem_nonlinear_fct(max_iter, eps, diff_const, S, Nx_elem, T, Nt_elem, ...
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
max_iter_TR = max_iter;
TR_rad_max = 1e9;
TR_rad = 0.75;
eta = 0.2;

for iter=1:1 %max_iter
    % entire looping happens in TR_step atm ... 

% computing gradient and hessian    
grad_J = grad_J_eval(hx, Nx_elem, ht, Nt_elem, bdy_cond, sigma, u, J_ss_lin, J_su_lin, J_us_lin, J_uu_lin, c1);
Hess_J = Hess_J_eval(hx, Nx_elem, ht, Nt_elem, bdy_cond, sigma, u, J_ss_lin, J_su_lin, J_us_lin, J_uu_lin, c1);

Hess_J_inner = Hess_J([inner_ind_sigma, tot_pts+inner_ind_u], [inner_ind_sigma, tot_pts+inner_ind_u]);
grad_J_inner = grad_J([inner_ind_sigma, tot_pts+inner_ind_u]); 

%%%%%%%%%%%%%%%%%%%%%%%% MULTIGRID %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
update = zeros(2*tot_pts,1);

update_inner = - Hess_J_inner \ grad_J_inner;

update([inner_ind_sigma, tot_pts+inner_ind_u])=update_inner;

% TRUST REGION STEP
[sigma, u, norm_grad_J, ~, TR_rad, conv_data] = TR_step(TR_rad, TR_rad_max, eta, max_iter_TR, S, Nx_elem, T, Nt_elem, u0, bdy_cond, bdy_left, bdy_right, sigma, u, update, grad_J, Hess_J, J_ss_lin, J_su_lin, J_us_lin, J_uu_lin, c1, c2, diff_const, iter);
 

file_name = ['conv_data_', num2str(S), 'by', num2str(T),'and' num2str(Nx_elem), 'by', num2str(Nt_elem), 'elem.mat'];
save(file_name, 'conv_data');
%     if(norm_grad_J < eps)
%         fprintf('convergence criterion reached after %d iterations\n', iter);
%         break;
%     end
 
end
% ************************* PLOTTING ****************************** % 

% reshape
u_mat = zeros(Nx_pts, Nt_pts);
sigma_mat = zeros(Nx_pts, Nt_pts);

for ts=1:Nt_pts
    u_mat(:, ts) = u((ts-1)*Nx_pts+1:ts*Nx_pts);
    sigma_mat(:, ts) = sigma((ts-1)*Nx_pts+1:ts*Nx_pts);

end

font_size = 15;

figure;
mesh(x_vec, t_vec, transpose(u_mat));
xlabel('space', 'interpreter','latex', 'FontSize', font_size);
ylabel('time', 'interpreter','latex', 'FontSize', font_size);
zlabel('u','interpreter','latex', 'FontSize', font_size);
colormap winter;

%zlim([min(u), max(u)]);
%title('LSFEM direct solve, f(u) = u*(1-u)*(u-0.1)');


% figure
% mesh(x_vec, t_vec(1:end), transpose(sigma_mat(:,1:end)));
% xlabel('space');
% ylabel('time');
% zlabel('sigma');
% %zlim([min(u), max(u)]);
% title('SIGMA heat equation LSFEM nonlinear');

end
