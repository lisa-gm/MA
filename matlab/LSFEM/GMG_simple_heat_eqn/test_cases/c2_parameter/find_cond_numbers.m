% ******************************************************** %
% ********** PLOTTING COMPARISON VALUES C2 *************** %
% ******************************************************** %

% ********** DIFFERENCES IN DIFFUSION CONSTANT *********** %

close all;
clear all;
% load solutions
path = '/Users/usi/Dropbox/MA/matlab/LSFEM/GMG_simple_heat_eqn/test_cases/c2_parameter/';

% NOW DIFFUSION CONSTANT D = 10^(-2)

% #1
load([path, 'data_20by20_fine_elem_c2_1_d_1_10_-2.mat']);
sol_direct_c2_1 = sol_direct;
sol_mg_c2_1 = it_sol_mg;

it_res_mg_c2_1 = it_res_mg;
conv_coeff_res_c2_1 = conv_coeff_res;

Hess_J_c2_1 = Hess_J_sym;
cond(Hess_J_c2_1)

% #10
load('data_20by20_fine_elem_c2_10_d_1_10_-2.mat');
sol_direct_c2_10 = sol_direct;
sol_mg_c2_10 = it_sol_mg;

it_res_mg_c2_10 = it_res_mg;
conv_coeff_res_c2_10 = conv_coeff_res;

Hess_J_c2_10 = Hess_J_sym;
cond(Hess_J_sym)

% #100
load('data_20by20_fine_elem_c2_100_d_1_10_-2.mat');
sol_direct_c2_100 = sol_direct;
sol_mg_c2_100 = it_sol_mg;

conv_coeff_res_c2_100 = conv_coeff_res;
it_res_mg_c2_100 = it_res_mg;

Hess_J_c2_100 = Hess_J_sym;
cond(Hess_J_sym)


fprintf('geometric mean value res conv coeff c2 = 1  : %d \n', geomean(conv_coeff_res_c2_1));
fprintf('geometric mean value res conv coeff c2 = 10 : %d \n', geomean(conv_coeff_res_c2_10));
fprintf('geometric mean value res conv coeff c2 = 100: %d \n', geomean(conv_coeff_res_c2_100));



% *************** PLOTTING **************************** %
font_size = 20;

% convergence rates 
figure;
semilogx(1:length(conv_coeff_res_c2_1),conv_coeff_res_c2_1, 'b', 'LineWidth', 2); hold on;
semilogx(1:length(conv_coeff_res_c2_10),conv_coeff_res_c2_10, 'r', 'LineWidth', 2); hold on;
semilogx(1:length(conv_coeff_res_c2_100),conv_coeff_res_c2_100, 'g', 'LineWidth', 2); 
lgnd = legend({'$c_2 = 1$', '$c_2 = 10$', '$c_2 = 100$'},'Location','southeast','interpreter','latex');
temp = [lgnd; lgnd.ItemText];
set(temp, 'FontSize', 15);
set(temp, 'interpreter','latex');
xlabel('Iterations[$k$]', 'interpreter','latex', 'FontSize', font_size);
ylabel('Convergence Coefficient', 'interpreter','latex', 'FontSize', font_size);
saveas(gcf,'fig1.eps','epsc')

% residual over iterations
figure;
semilogy(1:length(it_res_mg_c2_1),it_res_mg_c2_1, 'b', 'LineWidth', 2); hold on;
semilogy(1:length(it_res_mg_c2_10),it_res_mg_c2_10, 'r', 'LineWidth', 2); hold on;
semilogy(1:length(it_res_mg_c2_100),it_res_mg_c2_100, 'g', 'LineWidth', 2); 
lgnd = legend('c2 = 1', 'c2 = 10', 'c2 = 100','Location','northeast');
temp = [lgnd; lgnd.ItemText];
set(temp, 'FontSize', 15);
set(temp, 'interpreter','latex');
xlabel('iterations', 'interpreter','latex', 'FontSize', font_size);
ylabel('norm residual', 'interpreter','latex', 'FontSize', font_size);



