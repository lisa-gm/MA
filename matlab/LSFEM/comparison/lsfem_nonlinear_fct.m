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
x_vec = linspace(int_space(1), int_space(2), Nx_pts);

int_time = [0,T];
Nt_pts = Nt_elem+1;
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

damp_fact = 0;
for iter=1:max_iter 
   damp_fact = (1+2*damp_fact)/3;

   [sigma, u, norm_grad_J, damp_fact] = newton_step_lsfem(S, Nx_elem, T, Nt_elem, u0, bdy_cond, bdy_left, bdy_right, sigma, u, J_ss_lin, J_su_lin, J_us_lin, J_uu_lin, c1, c2, diff_const, eps, damp_fact, iter);

    if(norm_grad_J < eps)
        fprintf('convergence criterion reached after %d iterations', iter);
        break;
    end
end

% ************************* PLOTTING ****************************** % 

% reshape
u_mat = zeros(Nx_pts, Nt_pts);
sigma_mat = zeros(Nx_pts, Nt_pts);

for ts=1:Nt_pts
    u_mat(:, ts) = u((ts-1)*Nx_pts+1:ts*Nx_pts);
    sigma_mat(:, ts) = sigma((ts-1)*Nx_pts+1:ts*Nx_pts);

end

save('u_mat.mat', 'u_mat');

figure
mesh(x_vec, t_vec(1:end), transpose(u_mat(:,1:end)));
xlabel('space');
ylabel('time');
zlabel('u');
%zlim([min(u), max(u)]);
title('heat equation LSFEM nonlinear');

% 
% figure
% mesh(x_vec, t_vec(1:end), transpose(sigma_mat(:,1:end)));
% xlabel('space');
% ylabel('time');
% zlabel('sigma');
% %zlim([min(u), max(u)]);
% title('SIGMA heat equation LSFEM nonlinear');

end
