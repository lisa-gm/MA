% ********************************************************************** %
% ****************** LOOKING AT EIGENFUNCTIONS ************************* %
% ********************************************************************** %

clear all;
close all;

% *********** DOMAIN SET UP *************** 
S = 1;
Nx_elem = 20;

T = 20;
Nt_elem = 40;

% ... then 

int_space = [0,S];
hx = (int_space(2) - int_space(1))/Nx_elem;
Nx_pts = Nx_elem+1;
x_vec = linspace(int_space(1), int_space(2), Nx_pts);

int_time = [0, T];
ht = (int_time(2) - int_time(1))/Nt_elem;
Nt_pts = Nt_elem+1;
t_vec = linspace(int_time(1), int_time(2), Nt_pts);

tot_elem = Nx_elem*Nt_elem;
tot_pts = Nx_pts*Nt_pts;

hxht = ht*hx;

% ********************* load Hessian of J ***************************** %
load('Hess_J.mat');
load('u_mat.mat');
hess_J_tot_pts = size(Hess_J, 1)/2;

% check if sizes match 
if(hess_J_tot_pts ~= tot_pts)
    fprintf('sizes dont match! \n');
    return;
end

% ********************* choose patches to look at *********************** %

% size 
no_pts_space = 3;
no_pts_time = 3;

% where in space & time in terms, assume that then following pts in sp + t
% look at several patches, constant areas, wave front ... 

% for constant 1 area
pts_in_sp_1 = 4;
pts_in_t_1 = 20;

% for constant 0 area
pts_in_sp_0 = 16;
pts_in_t_0 = 6;

% on wavefront
pts_in_sp_wf_a = 8;
pts_in_t_wf_a = 12;

pts_in_sp_wf_b = 11;
pts_in_t_wf_b = 34;

% ********************************************************************* %
% get local Hessians
[ind_list_1, Hess_J_loc_1, ind_list_0, Hess_J_loc_0, ind_list_wf_a, Hess_J_loc_wf_a, ind_list_wf_b, Hess_J_loc_wf_b] = get_local_Hess(Nx_pts, tot_pts, no_pts_space, ...
                no_pts_time, pts_in_sp_1, pts_in_t_1, pts_in_sp_0, pts_in_t_0, pts_in_sp_wf_a, pts_in_t_wf_a,  pts_in_sp_wf_b, pts_in_t_wf_b, Hess_J);

            
% ***************** visualise points on sol u *************************** %

figure(100); hold on;
mesh(x_vec, t_vec, u_mat'); hold on;
title('points chosen for submatrices');
xlabel('space');
ylabel('time');
for i=1:no_pts_space
    for j=1:no_pts_time
        plot3(x_vec(pts_in_sp_1+i-1),t_vec(pts_in_t_1+j-1),u_mat(pts_in_sp_1+i-1,pts_in_t_1+j-1),'.r','markersize',10);
        plot3(x_vec(pts_in_sp_0+i-1),t_vec(pts_in_t_0+j-1),u_mat(pts_in_sp_0+i-1,pts_in_t_0+j-1),'.r','markersize',10);
        
        plot3(x_vec(pts_in_sp_wf_a+i-1),t_vec(pts_in_t_wf_a+j-1),u_mat(pts_in_sp_wf_a+i-1,pts_in_t_wf_a+j-1),'.r','markersize',10);
        plot3(x_vec(pts_in_sp_wf_b+i-1),t_vec(pts_in_t_wf_b+j-1),u_mat(pts_in_sp_wf_b+i-1,pts_in_t_wf_b+j-1),'.r','markersize',10);
        
        %mesh(x(i), t(j), u_mat(i,j));
        hold on; 
    end
end


% get eigenfunctions & sort them
[V_1, D_1] = eig(Hess_J_loc_1); [~, perm_1] = sort(diag(D_1));
D_1 = D_1(perm_1, perm_1); V_1 = V_1(:, perm_1);

[V_0, D_0] = eig(Hess_J_loc_0); [~, perm_0] = sort(diag(D_0));
D_0 = D_0(perm_0, perm_0); V_0 = V_0(:, perm_0);

[V_wf_a, D_wf_a] = eig(Hess_J_loc_wf_a); [~, perm_wf_a] = sort(diag(D_wf_a));
D_wf_a = D_wf_a(perm_wf_a, perm_wf_a); V_wf_a = V_wf_a(:, perm_wf_a);


[V_wf_b, D_wf_b] = eig(Hess_J_loc_wf_b); [~, perm_wf_b] = sort(diag(D_wf_b));
D_wf_b = D_wf_b(perm_wf_b, perm_wf_b); V_wf_b = V_wf_b(:, perm_wf_b);

% save in cell
V = {V_1, V_0, V_wf_a, V_wf_b};
D = {D_1, D_0, D_wf_a, D_wf_b};

% ********************************************************************* %
% ****************** PLOTTING EIGENFUNCTIONS ************************** %
% ********************************************************************* %

no_eig_vec = size(V_1,2);

size_h = ceil(sqrt(no_eig_vec));
size_v = ceil(no_eig_vec/size_h);

descr = {'Constant 1 Area', 'Constant 0 Area', 'Wavefront a', 'Wavefront b'};

for s=1:4
    
    figure(s);
    suptitle(['Local Eigenfunctions of ', descr{s}]);
    for i=1:size_h
        for j=1:size_v
            if((i-1)*size_v+j <= no_eig_vec)
            subplot(size_v, size_h, (i-1)*size_v+j)
            surf(reshape(V{s}(:,(i-1)*size_v+j),2*no_pts_space,no_pts_time));
            %zlabel('V(:,1)');
            zlim([-0.5,0.5]);
            str = ['EigVal ', num2str(D{s}((i-1)*size_v+j, (i-1)*size_v+j))];
            title(str);
            end
        end
    end
end
