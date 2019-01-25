% ******************************************************** %
% **************** PLOTTING EIGENVALUES* ***************** %
% ******************************************************** %

% comparison between (1) "standard" case (2) works very well 
% (3) doesn't work well at all

% all of the same grid size.

% choose sample 2, 5, 11, 20

close all;
clear all;
 
% load solutions
path = '/Users/usi/Dropbox/MA/matlab/LSFEM/GMG_simple_heat_eqn/making_plots/eigenvalues/';

% all 400 elements, that is 20 by 20 
% 16 elements total 

tot_pts = 400;

% 2 
load([path, 'eigenvalues_test_case_2']);
D_it_mat_tc_2 = D_it_mat;

% 5 
load([path, 'eigenvalues_test_case_5.mat']);
D_it_mat_tc_5 = D_it_mat;


% 11
load([path, 'eigenvalues_test_case_11.mat']);
D_it_mat_tc_11 = D_it_mat;


% 20
load([path, 'eigenvalues_test_case_20.mat']);
D_it_mat_tc_20 = D_it_mat;


% ********************** PLOTTING ************************ %


font_size = 20; 

figure;
scatter(1:length(D_it_mat_tc_2), D_it_mat_tc_2, '.', 'b'); hold on;
scatter(1:length(D_it_mat_tc_5), D_it_mat_tc_5, '.', 'g'); hold on;
scatter(1:length(D_it_mat_tc_11), D_it_mat_tc_11, '.', 'r'); hold on;
scatter(1:length(D_it_mat_tc_20), D_it_mat_tc_20, '.', 'k'); 
%set(gca,'yscale','log');
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