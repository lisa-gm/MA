% ******************************************************** %
% ********** PLOTTING COMPARISON VALUES C2 *************** %
% ******************************************************** %

close all;
 
% load solutions

load('solutions_c2_1.mat');
sol_direct_c2_1 = sol_direct;
sol_mg_c2_1 = sol_mg;

load('solutions_c2_10.mat');
sol_direct_c2_10 = sol_direct;
sol_mg_c2_10 = sol_mg;

load('solutions_c2_100.mat');
sol_direct_c2_100 = sol_direct;
sol_mg_c2_100 = sol_mg;

% load residuals

load('it_res_mg_20by20_fine_elem_c2_1.mat');
it_res_mg_c2_1 = it_res_mg;

load('it_res_mg_20by20_fine_elem_c2_10.mat');
it_res_mg_c2_10 = it_res_mg;

load('it_res_mg_20by20_fine_elem_c2_100.mat');
it_res_mg_c2_100 = it_res_mg;

% load conv coeff residual

load('conv_coeff_res_20by20_fine_elem_c2_1.mat');
conv_coeff_res_c2_1 = conv_coeff_res;

load('conv_coeff_res_20by20_fine_elem_c2_10.mat');
conv_coeff_res_c2_10 = conv_coeff_res;

load('conv_coeff_res_20by20_fine_elem_c2_100.mat');
conv_coeff_res_c2_100 = conv_coeff_res;

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



