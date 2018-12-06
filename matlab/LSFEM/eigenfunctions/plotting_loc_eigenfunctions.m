% ********************************************************************** %
% ****************** LOOKING AT EIGENFUNCTIONS ************************* %
% ********************************************************************** %

% get info from what we are looking at

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
hess_J_tot_pts = size(Hess_J, 1)/2;

% check if sizes match 
if(hess_J_tot_pts ~= tot_pts)
    fprintf('sizes dont match! \n');
    return;
end

% ********************* choose patch to look at *********************** %

% size 
no_pts_space = 2;
no_pts_time = 2;

% where in space & time in terms, assume that then following pts in sp + t
% look at several patches, constant areas, wave front ... 
pt_in_sp = 5;
pt_in_t = 5;

% ********************************************************************* %
% compute indices in matrix, for sigma & u part
ind_list_sigma = [];
for ts=1:no_pts_time
    ind_list_sigma = [ind_list_sigma, Nx_pts*(pt_in_t-1)+(ts-1)*Nx_pts+(pt_in_sp:pt_in_sp+no_pts_space-1)]; 
end

ind_list_u = ind_list_sigma + tot_pts;
ind_list = [ind_list_sigma, ind_list_u];

% extract matrix entries
Hess_J_loc = Hess_J(ind_list, ind_list);

% get eigenfunctions
[V, D] = eig(Hess_J_loc);

% ********************************************************************* %
% ****************** PLOTTING EIGENFUNCTIONS ************************** %
% ********************************************************************* %

no_eig_vec = size(V,2);

size_h = ceil(sqrt(no_eig_vec));
size_v = ceil(no_eig_vec/size_h);

figure;
for i=1:size_h
    for j=1:size_v
        if((i-1)*size_v+j <= no_eig_vec)
        subplot(size_v, size_h, (i-1)*size_v+j)
        surf(reshape(V(:,(i-1)*size_v+j),2*no_pts_space,no_pts_time));
        %zlabel('V(:,1)');
        zlim([-0.5,0.5]);
        end
suptitle('Local Eigenfunctions');
    end
end
