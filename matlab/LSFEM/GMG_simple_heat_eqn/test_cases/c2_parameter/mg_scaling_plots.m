% ******************************************************** %
% ************** MULTIGRID SCALING PLOTS ***************** %
% ******************************************************** %

close all;
clear all;
 
% load solutions
path = '/Users/usi/Dropbox/MA/matlab/LSFEM/GMG_simple_heat_eqn/test_cases/c2_parameter/';

% scaling results for D = 0.1, c_2 = 10

% 16 elements total 
load([path, 'data_4by4_fine_elem_c2_10_d_1_10_-1.mat']);
it_res_mg_N4 = it_res_mg;
conv_coeff_res_N4 = conv_coeff_res;

% 36 elements total 
load([path, 'data_6by6_fine_elem_c2_10_d_1_10_-1.mat']);
it_res_mg_N6 = it_res_mg;
conv_coeff_res_N6 = conv_coeff_res;

% 100 elements total 
load([path, 'data_10by10_fine_elem_c2_10_d_1_10_-1.mat']);
it_res_mg_N10 = it_res_mg;
conv_coeff_res_N10 = conv_coeff_res;
cond(Hess_J_sym)

% 400 elements total 
load([path, 'data_20by20_fine_elem_c2_10_d_1_10_-1.mat']);
it_res_mg_N20 = it_res_mg;
conv_coeff_res_N20 = conv_coeff_res;
cond(Hess_J_sym)

% 900 elements total 
load([path, 'data_30by30_fine_elem_c2_10_d_1_10_-1.mat']);
it_res_mg_N30 = it_res_mg;
conv_coeff_res_N30 = conv_coeff_res;
cond(Hess_J_sym)

% 1600 elements total 
load([path, 'data_40by40_fine_elem_c2_10_d_1_10_-1.mat']);
it_res_mg_N40 = it_res_mg;
conv_coeff_res_N40 = conv_coeff_res;

% 1600 elements total 
load([path, 'data_46by46_fine_elem_c2_10_d_1_10_-1.mat']);
it_res_mg_N45 = it_res_mg;
conv_coeff_res_N45 = conv_coeff_res;

% 1600 elements total 
load([path, 'data_50by50_fine_elem_c2_10_d_1_10_-1.mat']);
it_res_mg_N50 = it_res_mg;
conv_coeff_res_N50 = conv_coeff_res;


% ************** PLOTS ************************************** %
font_size = 20;

% residual over iterations
% figure;
% semilogy(1:length(it_res_mg_N4),it_res_mg_N4, 'k', 'LineWidth', 2); hold on;
% semilogy(1:length(it_res_mg_N6),it_res_mg_N6, 'r', 'LineWidth', 2); hold on;
% semilogy(1:length(it_res_mg_N10),it_res_mg_N10, 'g', 'LineWidth', 2); hold on;
% semilogy(1:length(it_res_mg_N20),it_res_mg_N20, 'm', 'LineWidth', 2); hold on;
% semilogy(1:length(it_res_mg_N30),it_res_mg_N30, 'b', 'LineWidth', 2); hold on;
% semilogy(1:length(it_res_mg_N40),it_res_mg_N40, 'c', 'LineWidth', 2); hold on;
% 
% lgnd = legend('16   elem', '36   elem', '100  elem', '400  elem', '900  elem', '1600 elem', 'Location','northeast');
% temp = [lgnd; lgnd.ItemText];
% set(temp, 'FontSize', 15);
% set(temp, 'interpreter','latex');
% xlabel('iterations[k]', 'interpreter','latex', 'FontSize', font_size);
% ylabel('norm residual', 'interpreter','latex', 'FontSize', font_size);
% 
% % convergence rate
% figure;
% plot(1:length(conv_coeff_res_N4),conv_coeff_res_N4, 'k', 'LineWidth', 2); hold on;
% plot(1:length(conv_coeff_res_N6),conv_coeff_res_N6, 'r', 'LineWidth', 2); hold on;
% plot(1:length(conv_coeff_res_N10),conv_coeff_res_N10, 'g', 'LineWidth', 2); hold on;
% plot(1:length(conv_coeff_res_N20),conv_coeff_res_N20, 'm', 'LineWidth', 2); hold on;
% plot(1:length(conv_coeff_res_N30),conv_coeff_res_N30, 'b', 'LineWidth', 2); hold on;
% plot(1:length(conv_coeff_res_N40),conv_coeff_res_N40, 'c', 'LineWidth', 2); hold on;
% lgnd = legend('16   elem', '36   elem', '100  elem', '400  elem', '900  elem', '1600 elem', 'Location','northeast');
% temp = [lgnd; lgnd.ItemText];
% set(temp, 'FontSize', 15);
% set(temp, 'interpreter','latex');
% xlabel('iterations[k]', 'interpreter','latex', 'FontSize', font_size);
% ylabel('$\mu_k$', 'interpreter','latex', 'FontSize', font_size);
