function plot_sol_and_coarse_grid(Nx_c, Nt_c, x_vec_f, t_vec_f, u_mat, sigma_mat)
tot_pts = length(x_vec_f)*length(t_vec_f);

hx_f = x_vec_f(2) - x_vec_f(1);
ht_f = t_vec_f(2) - t_vec_f(1);

hx_c = (x_vec_f(end) - x_vec_f(1))/Nx_c;
ht_c = (t_vec_f(end) - t_vec_f(1))/Nt_c;

k_x = hx_c/hx_f;
k_t = ht_c/ht_f;

if (floor(k_x) ~= k_x || floor(k_t) ~= k_t)
    frintf('invalid coarse and fine grid ratio!');
    return;
end

ind_list_c_x = 1:k_x:length(x_vec_f);
ind_list_c_t = 1:k_t:length(t_vec_f);


sigma_mat_c = sigma_mat(ind_list_c_x, ind_list_c_t);
u_mat_c = u_mat(ind_list_c_x, ind_list_c_t);


figure; hold on;
mesh(x_vec_f, t_vec_f, sigma_mat'); hold on;
title('error in SIGMA + points on coarse grid');
xlabel('space');
ylabel('time');
for i=1:Nx_c + 1
    for j=1:Nt_c + 1
        plot3(x_vec_f(ind_list_c_x(i)),t_vec_f(ind_list_c_t(j)),sigma_mat_c(i,j),'.r','markersize',10);
        hold on; 
    end
end

figure; hold on;
mesh(x_vec_f, t_vec_f, u_mat'); hold on;
title('U + points on coarse grid');
xlabel('space');
ylabel('time');
for i=1:Nx_c + 1
    for j=1:Nt_c + 1
        plot3(x_vec_f(ind_list_c_x(i)),t_vec_f(ind_list_c_t(j)),u_mat_c(i,j),'.r','markersize',10);
        hold on; 
    end
end


% additionally label coarse grid points in plot 

% figure;
% subplot(1,2,1);
% mesh(x_vec_f, t_vec_f, transpose(sigma_mat_direct-sigma_mat_mg)); 
% xlabel('space');
% ylabel('time');
% zlabel('diff sigma');
% %zlim([min(u), max(u)]);
% title('sigma direct - sigma mg');
% 
% subplot(1,2,2);
% mesh(x_vec_f, t_vec_f, transpose(u_mat_direct-u_mat_mg)); 
% xlabel('space');
% ylabel('time');
% zlabel('diff u');
% %zlim([min(u), max(u)]);
% title('u direct - u mg');

end