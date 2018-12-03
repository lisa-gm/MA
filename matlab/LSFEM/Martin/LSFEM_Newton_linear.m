%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% LSFEM Newton linear %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;

max_iter = 1;

%% set up domain 
int_space = [0,1];
Nx_elem = 5;
hx = (int_space(2) - int_space(1))/Nx_elem;
Nx_pts = Nx_elem+1;
x_vec = linspace(int_space(1), int_space(2), Nx_pts);

int_time = [0,1];
Nt_elem = 5;
ht = (int_time(2) - int_time(1))/Nt_elem;
Nt_pts = Nt_elem+1;
t_vec = linspace(int_time(1), int_time(2), Nt_pts);

tot_elem = Nx_elem*Nt_elem;
tot_pts = Nx_pts*Nt_pts;

%% exact solution

u_exact = zeros(Nx_pts, Nt_pts);
sigma_exact = zeros(Nx_pts, Nt_pts);


for ts=1:Nt_pts
    %u_exact(:,ts) = cos(pi*t(ts))*sin(pi*x); % x.^2*t(ts); 
    %u_exact(:, ts) = 5;
    %u_exact(:, ts) = sin(2*pi*t_vec(ts)); %*cos(pi*t_vec(ts));
    %u_exact(:, ts) = t_vec(ts);
    %u_exact(:,ts) = x_vec;
    u_exact(:,ts) = (x_vec).^2;
    
    %sigma_exact(:, ts) = 2*pi*cos(2*pi*x_vec);
    sigma_exact(:,ts) = 2*x_vec;
    %sigma_exact(:, ts) = pi*cos(pi*x_vec)*cos(pi*t_vec(ts));
end

%% initializing arrays

f_sigma_int = zeros(tot_pts,1);
f_u_int = zeros(tot_pts,1);

test_f_sigma_int = zeros(tot_pts,1);

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

%% boundary conditions
bdy_cond = 'Dirichlet';
%bdy_cond = 'Neumann';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(strcmp(bdy_cond, 'Dirichlet'))
    bdy_left = 0*ones(Nt_pts, 1);
    bdy_right = 1*ones(Nt_pts, 1);
    %bdy_left = t_vec;
    %bdy_right = t_vec;
end

% assume neumann bdy cond of the form \sigma*n_l = g_l, n_l = -1
% \sigma*n_r = g_r, n_r = 1
% enforce this in the matrix
if(strcmp(bdy_cond, 'Neumann'))
    g_l = zeros(Nt_pts,1);
    g_r = zeros(Nt_pts,1);
    
    for ts=1:Nt_pts
        g_l(ts) = 0;
        g_r(ts) = 0;
    end  
end

%% initial conditions
%u0 = sin(pi*x_vec);
u0 = u_exact(:,1);
%u0 = 0;

u = ones(tot_pts, 1);
u(1:Nx_pts) = u0;
sigma = 0.0*ones(tot_pts,1);

u = reshape(u_exact, [tot_pts, 1]);
sigma = reshape(sigma_exact, [tot_pts, 1]);


u(1:Nx_pts) = u0;
u(1:Nx_pts:end) = bdy_left;
u(Nx_pts:Nx_pts:end) = bdy_right;

sol = [sigma; u];

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

% define polynomial here
a0= @(x,t) -2; %- sin(pi*t); %pi^2*sin(pi*x)*cos(pi*t) - pi*sin(pi*x)*sin(pi*t); 
a1= @(x,t) 0; 
a2= @(x,t) 0; 
a3=0;

f = @(u,x,t) a0(x,t) + a1*u + a2*u^2 + a3*u^3;
df = @(u) a1 + 2*a2*u + 3*a3*u^2;
d2f = @(u) 2*a2 + 6*a3*u;

% approximate for integration, should be exact for polynomials of deg 1
d2f_dep_xt = @(coeff_vec,x,t, pt_i, pt_j) 2*a2(x,t) + 6*a3*dot(coeff_vec, basis_fct_eval(qd_pts(pt_i), qd_pts(pt_j)));
df_dep_xt = @(coeff_vec,x,t, pt_i, pt_j) a1(x,t) + 2*a2(x,t)*dot(coeff_vec, basis_fct_eval(qd_pts(pt_i), qd_pts(pt_j))) + 3*a3*(dot(coeff_vec, basis_fct_eval(qd_pts(pt_i), qd_pts(pt_j))))^2;
f_dep_xt = @(coeff_vec,x,t, pt_i, pt_j) a0(x,t) + a1(x,t)*dot(coeff_vec, basis_fct_eval(qd_pts(pt_i), qd_pts(pt_j))) + ...
            a2(x,t)*(dot(coeff_vec, basis_fct_eval(qd_pts(pt_i), qd_pts(pt_j))))^2 + a3*(dot(coeff_vec, basis_fct_eval(qd_pts(pt_i), qd_pts(pt_j))))^3;
        
%% loop linear part 

for elem_j = 1:Nt_elem                                      % iterate through all elements
    for elem_i = 1:Nx_elem
        
     curr_ind = get_elem_ind(elem_i, elem_j);               % get vector with indices of current elem
     curr_elem = pos2id(elem_i, elem_j);                    % get current element ID
     
     curr_elem_x = int_space(1) + (elem_i-1)*hx;            % get global coordinate, lower left corner
     curr_elem_t = int_time(1) + (elem_j-1)*ht;             % of current elem
     
     coeff_vec_sigma = sigma(curr_ind);
     coeff_vec_u = u(curr_ind);
     
        for pt_j=1:qd_deg
            qd_pt_t = curr_elem_t + qd_pts(pt_j)*ht;            % get global quad pt in y, scaled

            for pt_i=1:qd_deg             
                qd_pt_x = curr_elem_x + qd_pts(pt_i)*hx;        % get global quad pt in x, scaled
                
                bf_eval = basis_fct_eval(qd_pts(pt_i),qd_pts(pt_j));   % eval basis fct in quad pts       
                grad_tau_x_eval = grad_x_fct_eval(qd_pts(pt_i),qd_pts(pt_j));
                grad_v_t_eval = grad_t_fct_eval(qd_pts(pt_i),qd_pts(pt_j));
                
                d2f_local = d2f_dep_xt(coeff_vec_u,qd_pt_x,qd_pt_t, pt_i, pt_j); 
                df_local = df_dep_xt(coeff_vec_u,qd_pt_x,qd_pt_t, pt_i, pt_j);
                f_local = f_dep_xt(coeff_vec_u,qd_pt_x,qd_pt_t, pt_i, pt_j);
                
                sigma_x_local = dg_x(coeff_vec_sigma, qd_pts(pt_j));   % to compute sigma_x term
                u_t_local = dg_t(coeff_vec_u, qd_pts(pt_i));
                
                % iterate though basis functions for each element

                for bf_j=1:4
                    
                    % compute rhs integrals for f = f(x,t)
                    f_sigma_int(curr_ind(bf_j)) = f_sigma_int(curr_ind(bf_j)) + ht*hx*grad_tau_x_eval(bf_j)*f_local*qd_weights(pt_i)*qd_weights(pt_j);            %ass rhs                                                  %rhs
                    f_u_int(curr_ind(bf_j)) = f_u_int(curr_ind(bf_j)) - ht*hx*grad_v_t_eval(bf_j)*f_local*qd_weights(pt_i)*qd_weights(pt_j);
                    
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

%% Newton iteration 
max_iter =2;
solutions = zeros(max_iter, tot_pts);

for iter=1:max_iter
    
%%%%%%%%%%% putting things together for gradient J %%%%%%%%%%%%%%

J_s = J_ss_lin*sigma + J_su_lin*u + f_sigma_int;
J_u = J_us_lin*sigma + J_uu_lin*u + f_u_int;

components_J_s = [J_ss_lin*sigma, J_su_lin*u, f_sigma_int];
components_J_u = [J_us_lin*sigma, J_uu_lin*u, f_u_int];

J_u(1:Nx_pts) = 0;
J_u(1:Nx_pts:end) = 0;
J_u(Nx_pts:Nx_pts:end) = 0;

grad_J = [J_s; J_u];

fprintf('norm residual: %d, iteration: %d\n ', norm(grad_J), iter);

%%%%%%%%%%%%%%%%%% putting things together ----- HESSIAN --------
% linear here

J_ss = J_ss_lin;
J_uu = J_uu_lin;

J_us = J_us_lin;
J_su = J_su_lin;

%%%%%%%%%%%% putting in Dirichlet boundary conditions

J_uu(1:Nx_pts,:) = 0;
J_uu(1:Nx_pts:end) = 0;
J_uu(Nx_pts:Nx_pts:end) = 0;

J_us(1:Nx_pts,:) = 0;
J_us(1:Nx_pts:end) = 0;
J_us(Nx_pts:Nx_pts:end) = 0;

J_uu(1:Nx_pts, 1:Nx_pts) = eye(Nx_pts, Nx_pts);
J_uu(1:Nx_pts:end, 1:Nx_pts:end) = eye(Nt_pts, Nt_pts);
J_uu(Nx_pts:Nx_pts:end, Nx_pts:Nx_pts:end) = eye(Nt_pts, Nt_pts);

Hess_J = [J_ss, J_su; J_us, J_uu];
block_size = size(J_ss,1);

%%%%%%%% SOLVE 
sol = sol - Hess_J \ grad_J;

sigma = sol(1:block_size);
u = sol(block_size+1:end);

u(1:Nx_pts) = u0;
u(1:Nx_pts:end) = bdy_left;
u(Nx_pts:Nx_pts:end) = bdy_right;

sol = [sigma; u];

solutions(iter, :) = u;

end
%% PLOT 

% reshape
u_mat = zeros(Nx_pts, Nt_pts);
u_mat_1 = zeros(Nx_pts, Nt_pts);


for ts=1:Nt_pts
    u_mat(:, ts) = u((ts-1)*Nx_pts+1:ts*Nx_pts);
    u_mat_1(:, ts) = solutions(1,(ts-1)*Nx_pts+1:ts*Nx_pts);
end

figure
mesh(x_vec, t_vec(1:end), transpose(u_mat(:,1:end)));
xlabel('space');
ylabel('time');
zlabel('u');
zlim([min(u), max(u)]);
title('heat equation LSFEM linear');

norm(u_mat - u_exact);

figure
mesh(x_vec, t_vec(1:end), transpose(abs(u_mat_1(:,1:end)-u_exact(:,1:end))));
xlabel('space');
ylabel('time');
zlabel('u');
%zlim([min(u), max(u)]);
title('difference u 1st iter u exact');