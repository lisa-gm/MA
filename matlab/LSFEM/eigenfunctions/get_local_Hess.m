% ********************************************************************* %
% compute indices in matrix, for sigma & u part

function [ind_list_1, Hess_J_loc_1, ind_list_0, Hess_J_loc_0, ind_list_wf_a, Hess_J_loc_wf_a, ind_list_wf_b, Hess_J_loc_wf_b] ...
                = get_local_Hess(Nx_pts, tot_pts, no_pts_space, no_pts_time, pts_in_sp_1, pts_in_t_1, pts_in_sp_0, pts_in_t_0, ...
                                 pts_in_sp_wf_a, pts_in_t_wf_a,  pts_in_sp_wf_b, pts_in_t_wf_b, Hess_J)
                             
ind_list_sigma_1 = []; ind_list_sigma_0 = []; 
ind_list_sigma_wf_a = []; ind_list_sigma_wf_b = [];

for ts=1:no_pts_time
    ind_list_sigma_1 = [ind_list_sigma_1, Nx_pts*(pts_in_t_1-1)+(ts-1)*Nx_pts+(pts_in_sp_1:pts_in_sp_1+no_pts_space-1)]; 
    ind_list_sigma_0 = [ind_list_sigma_0, Nx_pts*(pts_in_t_0-1)+(ts-1)*Nx_pts+(pts_in_sp_0:pts_in_sp_0+no_pts_space-1)]; 
    ind_list_sigma_wf_a = [ind_list_sigma_wf_a, Nx_pts*(pts_in_t_wf_a-1)+(ts-1)*Nx_pts+(pts_in_sp_wf_a:pts_in_sp_wf_a+no_pts_space-1)]; 
    ind_list_sigma_wf_b = [ind_list_sigma_wf_b, Nx_pts*(pts_in_t_wf_b-1)+(ts-1)*Nx_pts+(pts_in_sp_wf_b:pts_in_sp_wf_b+no_pts_space-1)]; 

end

ind_list_u_1 = ind_list_sigma_1 + tot_pts;
ind_list_1 = [ind_list_sigma_1, ind_list_u_1];

ind_list_u_0 = ind_list_sigma_0 + tot_pts;
ind_list_0 = [ind_list_sigma_0 ind_list_u_0];

ind_list_u_wf_a = ind_list_sigma_wf_a + tot_pts;
ind_list_wf_a = [ind_list_sigma_wf_a, ind_list_u_wf_a];

ind_list_u_wf_b = ind_list_sigma_wf_b + tot_pts;
ind_list_wf_b = [ind_list_sigma_wf_b, ind_list_u_wf_b];

% extract matrix entries
Hess_J_loc_1 = Hess_J(ind_list_1, ind_list_1);
Hess_J_loc_0 = Hess_J(ind_list_0, ind_list_0);
Hess_J_loc_wf_a = Hess_J(ind_list_wf_a, ind_list_wf_a);
Hess_J_loc_wf_b = Hess_J(ind_list_wf_b, ind_list_wf_b);