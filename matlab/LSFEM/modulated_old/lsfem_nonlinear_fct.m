%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% LSFEM NONLINEAR RHS %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%% PARAMETERS REMOVED %%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% ADDITIONAL FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [sigma_mat, u_mat] = lsfem_nonlinear_fct(max_iter, omega, eps, alpha, diff_const, S, Nx_elem, T, Nt_elem, ...
    u_exact, sigma_exact, bdy_cond, bdy_left, bdy_right, u0, a0, a1, a2, a3, c1, c2)

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
u = 0*ones(tot_pts, 1);
u = reshape(u0'*ones(1,Nt_pts), [tot_pts,1]);

%[sigma, u] = get_init_guess_lsfem(S, Nx_elem, T, Nt_elem, u0, bdy_cond, bdy_left, bdy_right, a0, a1, J_ss_lin, J_su_lin, J_us_lin, J_uu_lin, c1, c2, diff_const);


%%%%%%%%%%%%%%%%%%%%%%%%% NEWTON ITERATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
solutions = zeros(max_iter, tot_pts);

for iter=1:max_iter 
    
   [sigma, u, norm_grad_J] = newton_step_lsfem(S, Nx_elem, T, Nt_elem, u0, bdy_cond, bdy_left, bdy_right, sigma, u, a0, a1, a2, a3, alpha, J_ss_lin, J_su_lin, J_us_lin, J_uu_lin, c1, c2, diff_const, eps, iter);
    
   u_mat = zeros(Nx_pts, Nt_pts);
sigma_mat = zeros(Nx_pts, Nt_pts);

for ts=1:Nt_pts
    u_mat(:, ts) = u((ts-1)*Nx_pts+1:ts*Nx_pts);
    sigma_mat(:, ts) = sigma((ts-1)*Nx_pts+1:ts*Nx_pts);

end

% figure
% mesh(x_vec, t_vec(1:end), transpose(u_mat(:,1:end)));
% xlabel('space');
% ylabel('time');
% zlabel('u');
% %zlim([min(u), max(u)]);
% title('heat equation LSFEM nonlinear');
% 
% figure
% mesh(x_vec, t_vec(1:end), transpose(sigma_mat(:,1:end)));
% xlabel('space');
% ylabel('time');
% zlabel('sigma');
% %zlim([min(u), max(u)]);
% title('SIGMA heat equation LSFEM nonlinear');

   if(norm_grad_J < eps)
       break;
   end

solutions(iter, :) = u;
end

% ************************* PLOTTING ****************************** % 


% exact solution
% figure
% mesh(x, t(1:end), transpose(u_exact(:,1:end)));
% xlabel('space');
% ylabel('time');
% zlabel('u');
% %zlim([min(u_exact), max(u_exact)]);
% title('exact solution');

% reshape
u_mat = zeros(Nx_pts, Nt_pts);
sigma_mat = zeros(Nx_pts, Nt_pts);

for ts=1:Nt_pts
    u_mat(:, ts) = u((ts-1)*Nx_pts+1:ts*Nx_pts);
    sigma_mat(:, ts) = sigma((ts-1)*Nx_pts+1:ts*Nx_pts);

end

figure
mesh(x_vec, t_vec(1:end), transpose(u_mat(:,1:end)));
xlabel('space');
ylabel('time');
zlabel('u');
%zlim([min(u), max(u)]);
title('heat equation LSFEM nonlinear');

figure
mesh(x_vec, t_vec(1:end), transpose(sigma_mat(:,1:end)));
xlabel('space');
ylabel('time');
zlabel('sigma');
%zlim([min(u), max(u)]);
title('SIGMA heat equation LSFEM nonlinear');

% figure
% mesh(x_vec, t_vec(1:end), transpose(abs(u_mat_1(:,1:end)-u_exact(:,1:end))));
% xlabel('space');
% ylabel('time');
% zlabel('u');
% %zlim([min(u), max(u)]);
% title('difference u 1st iter u exact');

end
