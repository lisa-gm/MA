% *************************************************************** %
% **************** MULTIGRID TEST PROBLEM *********************** %
% *************************************************************** %

clear all; 
close all;

% define MULTIGRID PARAMETERS

levels = 3;
smoother = 'GaussSeidel';
max_iter = 20;

% take final Hessian of the converged problem and final gradient to 
% reconstruct solution u 

% setting up the DOMAIN
S = 1;
Nx_elem = 20;
Nx_pts = Nx_elem + 1;
hx = S/Nx_elem;
x_vec = linspace(0,S,Nx_pts); 

T = 20;
Nt_elem = 40;
ht = T/Nt_elem;
Nt_pts = Nt_elem + 1;
t_vec = linspace(0,T, Nt_pts);

hxht = ht*hx;
tot_pts = Nx_pts*Nt_pts;

% set boundary conditions
bdy_cond = 'Neumann';

if(strcmp(bdy_cond, 'Dirichlet'))
    inner_ind_u = [];
    for ind=2:Nt_pts
        inner_ind_u = [inner_ind_u, (ind-1)*Nx_pts+2:ind*Nx_pts-1];
    end
    inner_ind_sigma = 1:tot_pts;
end

if(strcmp(bdy_cond, 'Neumann'))
    inner_ind_sigma = [];
    for ind=1:Nt_pts
        inner_ind_sigma = [inner_ind_sigma, (ind-1)*Nx_pts+2:ind*Nx_pts-1];
    end
    inner_ind_u = Nx_pts+1:tot_pts;
end

% LOADING THE VARIABLES 
load('final_Hess_J_inner.mat');
load('final_grad_J_inner.mat');
load('final_sol.mat');

% start with initial guess for error/update
update = zeros(2*tot_pts,1);
update_inner = 0.5*ones(length(inner_ind_u) + length(inner_ind_sigma), 1);

% check if positive definite, should give p = 0 iff pd
[~,p] = chol(Hess_J_inner);
p

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% CALL MUTLIGRID %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 update_inner = V_cycle(Hess_J_inner, -grad_J_inner, update_inner, levels, max_iter, smoother);
 update([inner_ind_sigma, inner_ind_u]) = update_inner;
 norm(update)
 
 sol = sol + update;
 
 update_sigma = update(1:tot_pts);
 update_u = update(tot_pts+1:end);
 
% PLOTTING 

% reshape
update_u_mat = zeros(Nx_pts, Nt_pts);

for ts=1:Nt_pts
    update_u_mat(:, ts) = update_u((ts-1)*Nx_pts+1:ts*Nx_pts);
end
 
figure;
mesh(x_vec, t_vec(1:end), transpose(update_u_mat(:,1:end)));
xlabel('space');
ylabel('time');
zlabel('update');
%zlim([min(u), max(u)]);
title('Visualise Update');
 
 





