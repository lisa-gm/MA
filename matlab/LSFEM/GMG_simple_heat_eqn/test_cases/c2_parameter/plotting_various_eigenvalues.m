% ******************************************************** %
% ***** PLOTTING EIGENVALUES* OF ITERATION MATRIX ******** %
% ******************************************************** %

% comparison between (1) "standard" case (2) works very well 
% (3) doesn't work well at all

% all of the same grid size.

% choose sample 2, 5, 11, 20

close all;
clear all;
 
% load solutions
path = '/Users/usi/Dropbox/MA/matlab/LSFEM/GMG_simple_heat_eqn/test_cases/c2_parameter/';

% all 400 elements, that is 20 by 20 
% 16 elements total 

tot_pts = 400;

% 2 
load([path, 'data_20by20_fine_elem_c2_1_d_1_10_-1.mat']);
length(it_res_mg)
Hess_J_sym_c2_1_d_0_1 = Hess_J_sym;
cond_c2_1_d_0_1 = cond(Hess_J_sym_c2_1_d_0_1)
cond_c2_1_d_0_1_ss = cond(Hess_J_sym_c2_1_d_0_1(1:tot_pts,1:tot_pts))
cond_c2_1_d_0_1_uu = cond(Hess_J_sym_c2_1_d_0_1(tot_pts+1:end,tot_pts+1:end))

% 5 
load([path, 'data_20by20_fine_elem_c2_1_d_1_10_0.mat']);
length(it_res_mg)
Hess_J_sym_c2_10_d_1 = Hess_J_sym;
cond_c2_10_d_1 = cond(Hess_J_sym_c2_10_d_1)
cond_c2_10_d_1_ss = cond(Hess_J_sym_c2_10_d_1(1:tot_pts,1:tot_pts))
cond_c2_10_d_1_uu = cond(Hess_J_sym_c2_10_d_1(tot_pts+1:end,tot_pts+1:end))

% 11
load([path, 'data_20by20_fine_elem_c2_10_d_1_10_-1.mat']);
length(it_res_mg)
Hess_J_sym_c2_10_d_0_1 = Hess_J_sym;
cond_c2_10_d_0_1 = cond(Hess_J_sym_c2_10_d_0_1)
cond_c2_10_d_0_1_ss = cond(Hess_J_sym_c2_10_d_0_1(1:tot_pts,1:tot_pts))
cond_c2_10_d_0_1_uu = cond(Hess_J_sym_c2_10_d_0_1(tot_pts+1:end,tot_pts+1:end))

% 20
load([path, 'data_20by20_fine_elem_c2_100_d_1_10_-1.mat']);
length(it_res_mg)
Hess_J_sym_c2_100_d_0_1 = Hess_J_sym;
cond_c2_100_d_0_1 = cond(Hess_J_sym_c2_100_d_0_1)
cond_c2_100_d_0_1_ss = cond(Hess_J_sym_c2_100_d_0_1(1:tot_pts,1:tot_pts))
cond_c2_100_d_0_1_uu = cond(Hess_J_sym_c2_100_d_0_1(tot_pts+1:end,tot_pts+1:end))


% *********** EIGENVALUES w/ BOUNDARY CONDITIONS ********* %

[V_sym_c2_1_d_0_1, D_sym_c2_1_d_0_1] = eig(Hess_J_sym_c2_1_d_0_1);
D_sym_c2_1_d_0_1 = sort(diag(D_sym_c2_1_d_0_1));

[V_sym_c2_10_d_1, D_sym_c2_10_d_1] = eig(Hess_J_sym_c2_10_d_1);
D_sym_c2_10_d_1 = sort(diag(D_sym_c2_10_d_1));

[V_sym_c2_10_d_0_1, D_sym_c2_10_d_0_1] = eig(Hess_J_sym_c2_10_d_0_1);
D_sym_c2_10_d_0_1 = sort(diag(D_sym_c2_10_d_0_1));

[V_sym_c2_100_d_1, D_sym_c2_100_d_0_1] = eig(Hess_J_sym_c2_100_d_0_1);
D_sym_c2_100_d_0_1 = sort(diag(D_sym_c2_100_d_0_1));

% ********************** PLOTTING ************************ %


font_size = 20; 

figure;
scatter(1:length(D_sym_c2_1_d_0_1), D_sym_c2_1_d_0_1, '.', 'b'); hold on;
scatter(1:length(D_sym_c2_10_d_1), D_sym_c2_10_d_1, '.', 'g'); hold on;
scatter(1:length(D_sym_c2_10_d_0_1), D_sym_c2_10_d_0_1, '.', 'r'); hold on;
scatter(1:length(D_sym_c2_100_d_0_1), D_sym_c2_100_d_0_1, '.', 'k'); 
set(gca,'yscale','log');
%set(gca,'xscale','log');
ylabel('eigenvalues', 'interpreter','latex', 'FontSize', font_size);
lgnd = legend('case 2', 'case 5', 'case 11', 'case 20', 'Location','southeast');
temp = [lgnd; lgnd.ItemText];
set(temp, 'FontSize', 15);
set(temp, 'interpreter','latex');

% figure;
% scatter(1:length(D_sym_c2_10_d_1), D_sym_c2_10_d_1, '.', 'b');
% set(gca,'yscale','log');
% ylabel('eigenvalues', 'interpreter','latex', 'FontSize', font_size);
% lgnd = legend('$c_2 = 10$, $D = 1$', 'Location','southeast');
% temp = [lgnd; lgnd.ItemText];
% set(temp, 'FontSize', 15);
% set(temp, 'interpreter','latex');
% 
% figure;
% scatter(1:length(D_sym_c2_10_d_0), D_sym_c2_10_d_0, '.', 'b');
% set(gca,'yscale','log');
% ylabel('eigenvalues', 'interpreter','latex', 'FontSize', font_size);
% lgnd = legend('$c_2 = 10$, $D = 0$', 'Location','southeast');
% temp = [lgnd; lgnd.ItemText];
% set(temp, 'FontSize', 15);
% set(temp, 'interpreter','latex');
% 
% figure;
% scatter(1:length(D_sym_c2_100_d_0_1), D_sym_c2_100_d_0_1, '.', 'b');
% set(gca,'yscale','log');
% ylabel('eigenvalues', 'interpreter','latex', 'FontSize', font_size);
% lgnd = legend('$c_2 = 100$, $D = 0$', 'Location','southeast');
% temp = [lgnd; lgnd.ItemText];
% set(temp, 'FontSize', 15);
% set(temp, 'interpreter','latex');