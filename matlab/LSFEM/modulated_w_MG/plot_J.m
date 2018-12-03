% ******************************************************************* %
% ******************* PLOT ENERGY LANDSCAPE ************************* %
% ******************************************************************* %

function J_eval_vec = plot_J(hx, Nx_elem, ht, Nt_elem, sigma, u, c1, c2, diff_const, iter)
Nx_pts = Nx_elem + 1;
Nt_pts = Nt_elem + 1;

x = 0:hx:Nx_elem*hx;
t = 0:ht:Nt_elem*ht;

[~, J_eval_vec] = J_eval_lsfem(hx, Nx_elem, ht, Nt_elem, sigma, u, c1, c2, diff_const);

J_eval_mat = reshape(J_eval_vec, [Nx_pts, Nt_pts]);

titl = ['Energy Landscape of J after ', num2str(iter), ' iter'];
% close all

figure;
mesh(x, t, J_eval_mat');
xlabel('space');
ylabel('time');
title(titl);

end