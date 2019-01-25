% ******************************************************** %
% ********* LINEARISED HESSIAN AT THE SOLUTION *********** % 
% ******************************************************** %

close all;
clear all;
% load Hessian, check if really positive definite

max_iter = 10;
sample_size = 100;
% want to see what exactly? 
% apply smoother and then see that damped super quickly on constant areas? 
% not damped along the front?! 


% contains the variables 'Hess_J', 'grad_J', 'S', 'Nx_elem', 'T', 'Nt_elem','sigma', 'u', 'u0', 'bdy_cond', 'bdy_left', 'bdy_right','diff_const')
load('Hess_J_20by40elem_diff0_001.mat');

Nx_pts = Nx_elem + 1;
Nt_pts = Nt_elem + 1;
tot_pts = Nx_pts*Nt_pts;

x_vec = linspace(0,S, Nx_pts);
t_vec = linspace(0,T, Nt_pts);

% now use some initial guess and see if it converges to the right solution
% if we use zero rhs, only enforce boundary conditions

[rhs_sym, Hess_J_sym] = apply_bdy_cond(Nx_pts, Nt_pts, Hess_J, grad_J, bdy_cond, u0, bdy_left, bdy_right, 1);

% how about looking at :
% (I + P*Hess_J) ....?!

% create vector over whole space-time domain that is zero everywhere 
% except for local patches of pertubations. 

p_size_s = 3;
p_size_t = 10;

% index lower left corner
index_s_const = 4;
index_t_const = 23;

index_s_wf =10;
index_t_wf = 20;

% now compute all other indices
ind_list_const = zeros(p_size_t, p_size_s);
ind_list_wf = zeros(p_size_t, p_size_s);


for it=1:p_size_t
    ind_list_const(it,:) = index_t_const*Nx_pts + (it - 1)*Nx_pts + [index_s_const:index_s_const+p_size_s-1];
    ind_list_wf(it,:) = index_t_wf*Nx_pts + (it - 1)*Nx_pts + [index_s_wf:index_s_wf+p_size_s-1];
end

% for now only pertubations in u?! what do we do with \sigma?
init_e_const = zeros(2*Nx_pts*Nt_pts,1);
init_e_wf = zeros(2*Nx_pts*Nt_pts,1);


relative_red_const = zeros(1, sample_size);
relative_red_wf = zeros(1, sample_size);

norm_pert_const = zeros(1, sample_size);
norm_pert_wf = zeros(1, sample_size);

grad_J = zeros(size(Hess_J,1),1);

norm_res_e_const = zeros(sample_size, max_iter+1);
norm_res_e_wf = zeros(sample_size, max_iter+1);

rho_list_const = zeros(1, sample_size);
c_list_const =  zeros(1, sample_size);

rho_list_wf = zeros(1, sample_size);
c_list_wf =  zeros(1, sample_size);

% save initial pertubations to make results reproduceable
list_init_pert_const = zeros(2*tot_pts, sample_size);
list_init_pert_wf = zeros(2*tot_pts, sample_size);

for k=1:sample_size
    
pert = 0.1*rand(p_size_s*p_size_t,1);

init_e_const(ind_list_const+tot_pts) = pert;
init_e_wf(ind_list_wf+tot_pts) = pert;

% also perturb in sigma by the same amount !!
%init_e_const(ind_list_const) = pert;
%init_e_wf(ind_list_wf) = pert;

init_e_const_u = init_e_const(tot_pts+1:end);
init_e_wf_u = init_e_wf(tot_pts+1:end);

list_init_pert_const(:, k) = init_e_const;
list_init_pert_wf(:, k) = init_e_wf;

[e_const, norm_res_e_const(k,:)] = JacobiSolve_LS(Hess_J, grad_J, init_e_const, max_iter);
[e_wf, norm_res_e_wf(k,:)] = JacobiSolve_LS(Hess_J, grad_J, init_e_wf, max_iter);

[c_list_const(k), rho_list_const(k)] = lin_regress(max_iter, norm_res_e_const(k,:));
[c_list_wf(k), rho_list_wf(k)] = lin_regress(max_iter, norm_res_e_wf(k,:));

u_e_const = e_const(Nx_pts*Nt_pts+1:end);
u_e_wf = e_wf(Nx_pts*Nt_pts+1:end);

sigma_e_const = e_const(1:tot_pts);
sigma_e_wf = e_wf(1:tot_pts);

% ***** HOW WELL DOES THE DAMPING WORK WHERE *********************** %
% looking at residuals

fprintf('initial residual const zone : %d\n', norm(grad_J - Hess_J*init_e_const));
fprintf('residual after %d iter const zone:  %d\n', max_iter, norm(grad_J - Hess_J*e_const));

fprintf('initial residual wf: %d\n', norm(grad_J - Hess_J*init_e_wf));
fprintf('residual after %d iter wf:  %d\n', max_iter, norm(grad_J - Hess_J*e_wf));

% compute relative reduction on the wave front and elsewhere 
relative_red_const(k) = norm(grad_J - Hess_J*(init_e_const))/ norm(grad_J - Hess_J*(e_const)) ;
relative_red_wf(k) = norm(grad_J - Hess_J*(init_e_wf))/ norm(grad_J - Hess_J*(e_wf));

end

save('list_init_pert_const.mat', 'list_init_pert_const');
save('list_init_pert_wf.mat', 'list_init_pert_wf');

fprintf('avg relative reduction of norm(res_init) / norm(res_final) constant zone : %d\n', mean(relative_red_const));
fprintf('avg relative reduction of  norm(res_init) / norm(res_final) wavefront : %d\n', mean(relative_red_wf));


% ************ DIRECT SOLVE TO CHECK IT WORKS *********************** %
% s_new = s_old - Hess_J(s_old)^_1 * grad_J(s_old)
% need s_old to obtain s_new 

sol_update = - Hess_J \ grad_J;

%sol_update = Hess_J_sym \ rhs_sym;

% compute sol = sol - Hess_J(sol) * grad_J(sol)
% since grad_J(sol) = 0;
% so what do i vary now? grad_J? to show that Hess_J^(-1) damps there quickly?

% look at smoother, ie
% (I - P*(Hess_J(sol)-grad_J)*e where e has pertubations in different spots ... 

% sol_new = sol_old - Hess_J^-1*grad_J
% --> update = sol_new - sol_old =  - Hess_J^-1*grad_J

% iteration would be
% Hess_J*sol_new = Hess_J*sol_old - grad_J --> but then we supp grad_J = 0


sigma_update = sol_update(1:Nx_pts*Nt_pts,1);
u_update = sol_update(Nx_pts*Nt_pts+1:end,1);

% ******************** PLOTTING SOLUTION **************************** %

% reshape
u_mat_pert_const = zeros(Nx_pts, Nt_pts);
u_mat_pert_init_const = zeros(Nx_pts, Nt_pts);

u_mat_pert_wf = zeros(Nx_pts, Nt_pts);
u_mat_pert_init_wf = zeros(Nx_pts, Nt_pts);

u_mat = zeros(Nx_pts, Nt_pts);
u_mat_update = zeros(Nx_pts, Nt_pts);

sigma_mat_pert_const = zeros(Nx_pts, Nt_pts);
sigma_mat_pert_wf = zeros(Nx_pts, Nt_pts);


for ts=1:Nt_pts
    %sigma_mat_direct(:,ts) = sigma_direct((ts-1)*Nx_pts+1:ts*Nx_pts);
    u_mat(:, ts) = u((ts-1)*Nx_pts+1:ts*Nx_pts);
    u_mat_update(:, ts) = u_update((ts-1)*Nx_pts+1:ts*Nx_pts);
    
    u_mat_pert_const(:, ts) = u_e_const((ts-1)*Nx_pts+1:ts*Nx_pts);
    u_mat_pert_init_const(:, ts) = init_e_const_u((ts-1)*Nx_pts+1:ts*Nx_pts);

    u_mat_pert_wf(:, ts) = u_e_wf((ts-1)*Nx_pts+1:ts*Nx_pts);
    u_mat_pert_init_wf(:, ts) = init_e_wf_u((ts-1)*Nx_pts+1:ts*Nx_pts);
    
    sigma_mat_pert_const(:, ts) = sigma_e_const((ts-1)*Nx_pts+1:ts*Nx_pts);
    sigma_mat_pert_wf(:, ts) = sigma_e_wf((ts-1)*Nx_pts+1:ts*Nx_pts);
end

font_size = 20;

figure;
subplot(2,2,1);
mesh(x_vec, t_vec, transpose(u_mat_pert_init_const));
xlabel('space', 'interpreter','latex', 'FontSize', font_size);
ylabel('time', 'interpreter','latex', 'FontSize', font_size);
zlabel('u', 'interpreter','latex', 'FontSize', font_size);
zlim([-0.02,0.1]);
title('initial pertubation constant zone', 'interpreter','latex', 'FontSize', font_size);

subplot(2,2,2);
mesh(x_vec, t_vec, transpose(u_mat_pert_const)); 
xlabel('space', 'interpreter','latex', 'FontSize', font_size);
ylabel('time', 'interpreter','latex', 'FontSize', font_size);
zlabel('u', 'interpreter','latex', 'FontSize', font_size);
zlim([-0.02,0.1]);
title('pertubation after smoothing constant zone', 'interpreter','latex', 'FontSize', font_size);

subplot(2,2,3);
mesh(x_vec, t_vec, transpose(u_mat_pert_init_wf));
xlabel('space', 'interpreter','latex', 'FontSize', font_size);
ylabel('time', 'interpreter','latex', 'FontSize', font_size);
zlabel('u', 'interpreter','latex', 'FontSize', font_size);
zlim([-0.02,0.1]);
title('initial pertubation wavefront', 'interpreter','latex', 'FontSize', font_size);

subplot(2,2,4);
mesh(x_vec, t_vec, transpose(u_mat_pert_wf)); 
xlabel('space', 'interpreter','latex', 'FontSize', font_size);
ylabel('time', 'interpreter','latex', 'FontSize', font_size);
zlabel('u', 'interpreter','latex', 'FontSize', font_size);
zlim([-0.02,0.1]);
title('pertubation after smoothing wavefront', 'interpreter','latex', 'FontSize', font_size);

figure;
plot(0:max_iter, norm_res_e_const(end,:), 'b', 0:max_iter, norm_res_e_wf(end,:), 'r', 'LineWidth',2); hold on;
lgnd = legend('on constant zone', 'on wavefront');
temp = [lgnd; lgnd.ItemText];
set(temp, 'FontSize', 15);
set(temp, 'interpreter','latex');
% title('residual norms pertubations');
xlabel('iterations', 'interpreter','latex', 'FontSize', font_size);
ylabel('residual norm', 'interpreter','latex', 'FontSize', font_size);


figure;
mesh(x_vec, t_vec, transpose(u_mat)); hold on;
for i=1:p_size_s
    for j=1:p_size_t
        plot3(x_vec(index_s_const+i-1),t_vec(index_t_const+j-1),u_mat(index_s_const+i-1,index_t_const+j-1),'.b','markersize',20);
        plot3(x_vec(index_s_wf+i-1),t_vec(index_t_wf+j-1),u_mat(index_s_wf+i-1,index_t_wf+j-1),'.r','markersize',20);
        
        %plot3(x_vec(pts_in_sp_wf_a+i-1),t_vec(pts_in_t_wf_a+j-1),u_mat(pts_in_sp_wf_a+i-1,pts_in_t_wf_a+j-1),'.r','markersize',10);
        %plot3(x_vec(pts_in_sp_wf_b+i-1),t_vec(pts_in_t_wf_b+j-1),u_mat(pts_in_sp_wf_b+i-1,pts_in_t_wf_b+j-1),'.r','markersize',10);
        
        %mesh(x(i), t(j), u_mat(i,j));
        hold on; 
    end
end
xlabel('space', 'interpreter','latex', 'FontSize', font_size);
ylabel('time', 'interpreter','latex', 'FontSize', font_size);
zlabel('u', 'interpreter','latex', 'FontSize', font_size);
colormap winter;
%zlim([-0.02,0.1]);


% and plot norms
fprintf('on avg const zone : | e^k | = %d * %d ^k\n', mean(c_list_const), mean(rho_list_const));
fprintf('on avg wf         : | e^k | = %d * %d ^k\n', mean(c_list_wf), mean(rho_list_wf));
