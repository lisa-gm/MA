%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% EVALUATE NONLINEAR TERMS %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TODO: deal with scaling from rest ...

% test with u = u(x), nonlinear f

clear all;
close all;

%% set up domain 
int_space = [0,1];
Nx_elem = 15;
hx = (int_space(2) - int_space(1))/Nx_elem;
Nx_pts = Nx_elem+1;
x = linspace(int_space(1), int_space(2), Nx_pts);

int_time = [0,1];
Nt_elem = 18;
ht = (int_time(2) - int_time(1))/Nt_elem;
Nt_pts = Nt_elem+1;
t = linspace(int_time(1), int_time(2), Nt_pts);

tot_elem = Nx_elem*Nt_elem;
tot_pts = Nx_pts*Nt_pts;

%% initializing arrays

f_sigma_int = zeros(tot_pts,1);
f_u_int = zeros(tot_pts,1);

% linear parts of Hessian
J_ss_lin = zeros(tot_pts);
J_su_lin = zeros(tot_pts);                          % compute both to check symmetry
J_us_lin = zeros(tot_pts);
J_uu_lin = zeros(tot_pts);

% for non linear parts of J_us, J_suX
tau_x_df_int = zeros(tot_pts);

% for non linear parts of J_uu
df_df_int = zeros(tot_pts);
f_d2f_int = zeros(tot_pts);
u_t_d2f_int = zeros(tot_pts);
s_x_d2f_int = zeros(tot_pts);

%% initial conditions
%u0 = - x.^2 + 2.*x + 1;
%u0 = sin(pi.*x) + 1;
u0 = 0;

%% boundary conditions
% bdy_left = ones(Nt_pts, 1);
% bdy_right = ones(Nt_pts, 1);

bdy_left = 6*t;
bdy_right = 6*t;

%% right hand side f
g = @(x,t) 6;

%% quadrature
qd_deg = 3;
qd_pts = 0.5*([-sqrt(3/5),0,+sqrt(3/5)]+1);                         %3-point gauss on [0,1]
qd_weights = 0.5*[5/9,8/9,5/9];

%% 

pos2id = @(i,j) (j-1)*Nx_pts + (i-1) + 1; 
get_elem_ind = @(i,j) [pos2id(i,j); pos2id(i+1,j); pos2id(i,j+1); pos2id(i+1,j+1)];

% to obtain sigma_x or u_x for iterate x_k = [sigma_k, u_k]
dg_x = @(coeff_vec, t) 1/hx*(coeff_vec(1)*(-1)*(1-t) + coeff_vec(2)*(1-t) + coeff_vec(3)*t*(-1) + coeff_vec(4)*t);
dg_t = @(coeff_vec, x) 1/ht*(coeff_vec(1)*(1-x)*(-1) + coeff_vec(2)*x*(-1) + coeff_vec(3)*(1-x) + coeff_vec(4)*x*(-1));

%% basis functions

basis_fct_eval = @(x,t) [(1-x)*(1-t); x*(1-t); (1-x)*t;  x*t];
grad_x_fct_eval = @(x,t) 1/hx*[(-1)*(1-t); 1*(1-t); (-1)*t;  1*t];
grad_t_fct_eval = @(x,t) 1/ht*[(1-x)*(-1); x*(-1); (1-x)*1;  x*1];

% function F=basis_function(x,t)
% F = [(1-x)*(1-t); x*(1-t); (1-x)*t;  x*t];
% end

%% define f
alpha = 0.1;

% define polynomial here
a0=0; a1=0; a2=0; a3=1/6;

f = @(u) a0 + a1*u + a2*u^2 + a3*u^3;
df = @(u) a1 + 2*a2*u + 3*a3*u^2;
d2f = @(u) 2*a2 + 6*a3*u;

% approximate for integration, should be exact for polynomials of deg 1
d2f_dep_xt = @(coeff_vec,x,t) 2*a2 + 6*a3*dot(coeff_vec, basis_fct_eval(x,t));
df_dep_xt = @(coeff_vec,x,t) a1 + 2*a2*dot(coeff_vec, basis_fct_eval(x,t)) + 6*a3*(dot(coeff_vec, basis_fct_eval(x,t)))^2;
f_dep_xt = @(coeff_vec,x,t) a0 + a1*dot(coeff_vec, basis_fct_eval(x,t)) + ...
            2*a2*(dot(coeff_vec, basis_fct_eval(x,t)))^2 + 6*a3*(dot(coeff_vec, basis_fct_eval(x,t)))^3;

%% assume that we have a current solution u_k, is a matrix, with values in each point
u_vec = ones(tot_pts,1);
sigma_vec = ones(tot_pts, 1);

%% loop

for elem_j = 1:Nt_elem                                      % iterate through all elements
    for elem_i = 1:Nx_elem
        
     curr_ind = get_elem_ind(elem_i, elem_j);               % get vector with indices of current elem
     curr_elem = pos2id(elem_i, elem_j);                    % get current element ID
     
     curr_elem_x = int_space(1) + (elem_i-1)*hx;            % get global coordinate, lower left corner
     curr_elem_t = int_time(1) + (elem_j-1)*ht;             % of current elem
     
     coeff_vec_sigma = sigma_vec(curr_ind);
     coeff_vec_u = u_vec(curr_ind);
     
        for pt_j=1:qd_deg
            qd_pt_t = curr_elem_t + qd_pts(pt_j)*ht;            % get global quad pt in y, scaled

            for pt_i=1:qd_deg             
                qd_pt_x = curr_elem_x + qd_pts(pt_i)*hx;        % get global quad pt in x, scaled
                
                bf_eval = basis_fct_eval(qd_pts(pt_i),qd_pts(pt_j));   % eval basis fct in quad pts       
                grad_tau_x_eval = grad_x_fct_eval(qd_pts(pt_i),qd_pts(pt_j));
                grad_v_t_eval = grad_t_fct_eval(qd_pts(pt_i),qd_pts(pt_j));
                
                d2f_local = d2f_dep_xt(coeff_vec_u,qd_pt_x,qd_pt_t); 
                df_local = df_dep_xt(coeff_vec_u,qd_pt_x,qd_pt_t);
                f_local = f_dep_xt(coeff_vec_u,qd_pt_x,qd_pt_t);
                
                g_local = g(qd_pt_x, qd_pt_t);
                
                sigma_x_local = dg_x(coeff_vec_sigma, qd_pt_t);        % to compute sigma_x term
                u_t_local = dg_t(coeff_vec_u, qd_pt_x);
                
                % iterate though basis functions for each element

                for bf_j=1:4
                    
                    % compute rhs integrals for f = f(x,t)
                    f_sigma_int(curr_ind(bf_j)) = f_sigma_int(curr_ind(bf_j)) - ht*hx*grad_tau_x_eval(bf_j)*g_local*qd_weights(pt_i)*qd_weights(pt_j);            %ass rhs                                                  %rhs
                    f_u_int(curr_ind(bf_j)) = f_u_int(curr_ind(bf_j)) + ht*hx*grad_v_t_eval(bf_j)*g_local*qd_weights(pt_i)*qd_weights(pt_j);
                    
                    for bf_i=1:4
                       % J_ss linear part: <\tau, \tau> + <\tau_x, \tau_x> 
                       J_ss_lin(curr_ind(bf_i), curr_ind(bf_j)) = J_ss_lin(curr_ind(bf_i), curr_ind(bf_j))+ht*hx*bf_eval(bf_i)*bf_eval(bf_j)*qd_weights(pt_i)*qd_weights(pt_j);
                       J_ss_lin(curr_ind(bf_i), curr_ind(bf_j)) = J_ss_lin(curr_ind(bf_i), curr_ind(bf_j))+ht*hx*grad_tau_x_eval(bf_i)*grad_tau_x_eval(bf_j)*qd_weights(pt_i)*qd_weights(pt_j);
                       
                       % J_su linear part: -<\tau, v_x> - <v_t, \tau_x
                       J_su_lin(curr_ind(bf_i), curr_ind(bf_j)) = J_su_lin(curr_ind(bf_i), curr_ind(bf_j)) - ht*hx*bf_eval(bf_i)*grad_tau_x_eval(bf_j)*qd_weights(pt_i)*qd_weights(pt_j);
                       J_su_lin(curr_ind(bf_i), curr_ind(bf_j)) = J_su_lin(curr_ind(bf_i), curr_ind(bf_j)) - ht*hx*grad_v_t_eval(bf_j)*grad_tau_x_eval(bf_i)*qd_weights(pt_i)*qd_weights(pt_j);
                       
                       % J_us linear part
                       J_us_lin(curr_ind(bf_i), curr_ind(bf_j)) = J_us_lin(curr_ind(bf_i), curr_ind(bf_j)) - ht*hx*bf_eval(bf_j)*grad_tau_x_eval(bf_i)*qd_weights(pt_i)*qd_weights(pt_j);
                       J_us_lin(curr_ind(bf_i), curr_ind(bf_j)) = J_us_lin(curr_ind(bf_i), curr_ind(bf_j)) - ht*hx*grad_v_t_eval(bf_i)*grad_tau_x_eval(bf_j)*qd_weights(pt_i)*qd_weights(pt_j);

                       % J_uu linear part:
                       J_uu_lin(curr_ind(bf_i), curr_ind(bf_j)) = J_uu_lin(curr_ind(bf_i), curr_ind(bf_j)) + ht*hx*grad_v_t_eval(bf_i)*grad_v_t_eval(bf_j)*qd_weights(pt_i)*qd_weights(pt_j);
                       J_uu_lin(curr_ind(bf_i), curr_ind(bf_j)) = J_uu_lin(curr_ind(bf_i), curr_ind(bf_j)) + ht*hx*grad_tau_x_eval(bf_i)*grad_tau_x_eval(bf_j)*qd_weights(pt_i)*qd_weights(pt_j);

                    end
                end
                   
            end
        end
     end
    
end

%% putting things together
J_ss = J_ss_lin;

%J_su = J_su_lin;  % + nonlinear bit later
J_uu = J_uu_lin;

% for now
J_us = J_us_lin;
J_su = J_su_lin;

%% putting in Dirichlet boundary conditions

J_uu(1:Nx_pts,:) = 0;
J_uu(1:Nx_pts:end) = 0;
J_uu(Nx_pts:Nx_pts:end) = 0;

J_us(1:Nx_pts,:) = 0;
J_us(1:Nx_pts:end) = 0;
J_us(Nx_pts:Nx_pts:end) = 0;

J_uu(1:Nx_pts, 1:Nx_pts) = eye(Nx_pts, Nx_pts);
J_uu(1:Nx_pts:end, 1:Nx_pts:end) = eye(Nt_pts, Nt_pts);
J_uu(Nx_pts:Nx_pts:end, Nx_pts:Nx_pts:end) = eye(Nt_pts, Nt_pts);

f_u_int(1:Nx_pts) = u0;
f_u_int(1:Nx_pts:end) = bdy_left;
f_u_int(Nx_pts:Nx_pts:end) = bdy_right;

% how to loop through several elements
% get indices from global vector

Hess_J = [J_ss, J_su; J_us, J_uu];
block_size = size(J_ss,1);

%% SOLVE 
f_vec_sigma = zeros(tot_pts,1);
F = [f_sigma_int; f_u_int];

sol = Hess_J \ F;
u = sol(block_size+1:end);

save('sol.mat', 'sol');
%% PLOT 

% reshape
u_mat = zeros(Nx_pts, Nt_pts);

for ts=1:Nt_pts
    u_mat(:, ts) = u((ts-1)*Nx_pts+1:ts*Nx_pts);
end

figure
mesh(x, t(1:end), transpose(u_mat(:,1:end)));
xlabel('space');
ylabel('time');
zlabel('u');
zlim([min(u), max(u)]);
title('heat equation LSFEM');