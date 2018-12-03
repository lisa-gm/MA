%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% LSFEM NONLINEAR RHS %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PUT IN SCALING FOR DIFFUSION, call diff_const_m

clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% GENERAL PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

max_iter = 5;                          % number of newton iterations
omega = 1;                              % for a scaled newton step
eps = 10^-12;

alpha = 1;                              % scaling in front of f(u) ie mu*f(u) -> derivatives
diff_const = 1;                         % diffusion constant, define through sigma = diff_const*nabla(u) 
                                        % hence doesn't appear in f terms..

c1 = 1;                                 % constants in front of J = c1*|| ....-f(u) ||^2 + c2*|| ...||^2
c2 = 1;                                 % need to multiply with them in computation of each integral 
                                        % so that one obtains correct
                                        % weighting for each one
                                        
%% DOMAIN SET UP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
int_space = [0,1];
Nx_elem = 20;
hx = (int_space(2) - int_space(1))/Nx_elem;
Nx_pts = Nx_elem+1;
x = linspace(int_space(1), int_space(2), Nx_pts);

int_time = [0,2];
Nt_elem = 40;
ht = (int_time(2) - int_time(1))/Nt_elem;
Nt_pts = Nt_elem+1;
t = linspace(int_time(1), int_time(2), Nt_pts);

tot_elem = Nx_elem*Nt_elem;
tot_pts = Nx_pts*Nt_pts;

hxht = ht*hx;
%% COMPUTING EXACT SOLUTION (if possible) %%%%%%%%%%%%%%%%%%
u_exact = zeros(Nx_pts, Nt_pts);

for ts=1:Nt_pts
    %u_exact(:,ts) = x.^2*t(ts); % cos(pi*t(ts))*sin(pi*x)
    u_exact(:, ts) = t(ts);
end

%% BOUNDARY CONDITIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%bdy_left = zeros(Nt_pts, 1);
%bdy_right = ones(Nt_pts, 1);
bdy_left = t;
bdy_right = t;

%% initial conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%u0 = -x.*(x-1)+1;
%u0 = sin(pi*x);
u0 = 0;

% u = ones(tot_pts, 1);
% 
% u(1:Nx_pts) = u0;
% u(1:Nx_pts:end) = bdy_left;
% u(Nx_pts:Nx_pts:end) = bdy_right;
% 
% sigma = ones(tot_pts,1);
% 
% sol = [sigma; u];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% ADDITIONAL FUNCTIONS %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% quadrature
qd_deg = 3;
qd_pts = 0.5*([-sqrt(3/5),0,+sqrt(3/5)]+1);                         %3-point gauss on [0,1]
qd_weights = 0.5*[5/9,8/9,5/9];

%% global indexing 

pos2id = @(i,j) (j-1)*Nx_pts + (i-1) + 1;                                           %to go from matrix to vector ID
get_elem_ind = @(i,j) [pos2id(i,j); pos2id(i+1,j); pos2id(i,j+1); pos2id(i+1,j+1)]; %give lower left matrix corner ID, get 
                                                                                    %vector IDs of the entire element

%% directional derivatives in direction x and t

dg_x = @(coeff_vec, t) 1/hx*(coeff_vec(1)*(-1)*(1-t) + coeff_vec(2)*(1-t) + coeff_vec(3)*t*(-1) + coeff_vec(4)*t);
dg_t = @(coeff_vec, x) 1/ht*(coeff_vec(1)*(1-x)*(-1) + coeff_vec(2)*x*(-1) + coeff_vec(3)*(1-x) + coeff_vec(4)*x);

%% basis functions, and their gradients

basis_fct_eval = @(x,t) [(1-x)*(1-t); x*(1-t); (1-x)*t;  x*t];
grad_x_fct_eval = @(x,t) 1/hx*[(-1)*(1-t); 1*(1-t); (-1)*t;  1*t];
grad_t_fct_eval = @(x,t) 1/ht*[(1-x)*(-1); x*(-1); (1-x)*1;  x*1];

% function F=basis_function(x,t)
% F = [(1-x)*(1-t); x*(1-t); (1-x)*t;  x*t];
% end

%% define f, polynomial of the form: sum_i a_i u^i with a_i = a(x,t), i<=3

%a0 = @(x,t) - pi*sin(pi*x)*sin(pi*t) + pi^2*sin(pi*x)*cos(pi*t) - (sin(pi*x)*cos(pi*t))^3;
a0 = @(x,t) 1 - t;
a1 = @(x,t) 1; 
a2 = @(x,t) 0; 
a3 = 0;
% 
f = @(u, x, t) alpha*(a0(x,t) + a1(x,t)*u + a2(x,t)*u^2 + a3*u^3);
df = @(u) alpha*(a1(x,t) + 2*a2(x,t)*u + 3*a3*u^2);
d2f = @(u) alpha*(2*a2(x,t) + 6*a3*u);

% approximate for integration, should be exact for polynomials of deg 1
d2f_dep_xt = @(coeff_vec,x,t) 2*a2(x,t) + 6*a3*dot(coeff_vec, basis_fct_eval(x,t));
df_dep_xt = @(coeff_vec,x,t) a1(x,t) + 2*a2(x,t)*dot(coeff_vec, basis_fct_eval(x,t)) + 3*a3*(dot(coeff_vec, basis_fct_eval(x,t)))^2;
f_dep_xt = @(coeff_vec,x,t) a0(x,t) + a1(x,t)*dot(coeff_vec, basis_fct_eval(x,t)) + ...
            a2(x,t)*(dot(coeff_vec, basis_fct_eval(x,t)))^2 + a3*(dot(coeff_vec, basis_fct_eval(x,t)))^3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% CONSRUCTING THE OPERATORS %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
%% ARRAY ALLOCATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% linear parts of Hessian
J_ss_lin = zeros(tot_pts);
J_su_lin = zeros(tot_pts);                          % compute both to check symmetry
J_us_lin = zeros(tot_pts);
J_uu_lin = zeros(tot_pts);

% additional linear terms for better initial guess
% solving ||u_t - div(\sigma) - a0 - a1*u ||^2 + || \sigma - grad(u) ||^2 first
% have additional terms including < ..., a0 + a1*u >

J_s_a0 = zeros(tot_pts,1);
J_u_a0 = zeros(tot_pts,1);

J_su_a1 = zeros(tot_pts);
J_us_a1 = zeros(tot_pts);
J_uu_a1 = zeros(tot_pts);
        
%% loop linear part, construct matrices

for elem_j = 1:Nt_elem                                      % iterate through all elements
    for elem_i = 1:Nx_elem
        
     curr_ind = get_elem_ind(elem_i, elem_j);               % get vector with indices of current elem
     curr_elem = pos2id(elem_i, elem_j);                    % get current element ID
     
     curr_elem_x = int_space(1) + (elem_i-1)*hx;            % get global coordinate, lower left corner
     curr_elem_t = int_time(1) + (elem_j-1)*ht;             % of current elem
     
        for pt_j=1:qd_deg                                   % iterate through quadrature points 
            qd_pt_t = curr_elem_t + qd_pts(pt_j)*ht;            % get global quad pt in y, scaled

            for pt_i=1:qd_deg             
                qd_pt_x = curr_elem_x + qd_pts(pt_i)*hx;        % get global quad pt in x, scaled
                
                bf_eval = basis_fct_eval(qd_pts(pt_i),qd_pts(pt_j));   % eval basis fct in quad pts       
                grad_tau_x_eval = grad_x_fct_eval(qd_pts(pt_i),qd_pts(pt_j));
                grad_v_t_eval = grad_t_fct_eval(qd_pts(pt_i),qd_pts(pt_j));
                
                                                      
                for bf_j=1:4                                           % iterate though basis functions
                              
                    % compute rhs integrals for f = f(x,t)
                    J_s_a0(curr_ind(bf_j)) = J_s_a0(curr_ind(bf_j)) - c1*hxht*grad_tau_x_eval(bf_j)*a0(qd_pt_x, qd_pt_t)*qd_weights(pt_i)*qd_weights(pt_j);            %ass rhs                                                  %rhs
                    J_u_a0(curr_ind(bf_j)) = J_u_a0(curr_ind(bf_j)) + c1*hxht*grad_v_t_eval(bf_j)*a0(qd_pt_x, qd_pt_t)*qd_weights(pt_i)*qd_weights(pt_j);
                    J_u_a0(curr_ind(bf_j)) = J_u_a0(curr_ind(bf_j)) - c1*hxht*bf_eval(bf_j)*a1(qd_pt_x, qd_pt_t)*a0(qd_pt_x, qd_pt_t)*qd_weights(pt_i)*qd_weights(pt_j);

                    
                    for bf_i=1:4
                       % J_ss linear part: <\tau, \tau> + <\tau_x, \tau_x> 
                       J_ss_lin(curr_ind(bf_i), curr_ind(bf_j)) = J_ss_lin(curr_ind(bf_i), curr_ind(bf_j)) + c2*hxht*bf_eval(bf_i)*bf_eval(bf_j)*qd_weights(pt_i)*qd_weights(pt_j);
                       J_ss_lin(curr_ind(bf_i), curr_ind(bf_j)) = J_ss_lin(curr_ind(bf_i), curr_ind(bf_j)) + c1*hxht*grad_tau_x_eval(bf_i)*grad_tau_x_eval(bf_j)*qd_weights(pt_i)*qd_weights(pt_j);
                       
                       % J_su linear part: -< \tau, v_x> - <v_t, \tau_x
                       J_su_lin(curr_ind(bf_i), curr_ind(bf_j)) = J_su_lin(curr_ind(bf_i), curr_ind(bf_j)) - c2*hxht*diff_const*bf_eval(bf_i)*grad_tau_x_eval(bf_j)*qd_weights(pt_i)*qd_weights(pt_j);
                       J_su_lin(curr_ind(bf_i), curr_ind(bf_j)) = J_su_lin(curr_ind(bf_i), curr_ind(bf_j)) - c1*hxht*grad_v_t_eval(bf_j)*grad_tau_x_eval(bf_i)*qd_weights(pt_i)*qd_weights(pt_j);
                       
                       % J_us linear part
                       J_us_lin(curr_ind(bf_i), curr_ind(bf_j)) = J_us_lin(curr_ind(bf_i), curr_ind(bf_j)) - c2*hxht*diff_const*bf_eval(bf_j)*grad_tau_x_eval(bf_i)*qd_weights(pt_i)*qd_weights(pt_j);
                       J_us_lin(curr_ind(bf_i), curr_ind(bf_j)) = J_us_lin(curr_ind(bf_i), curr_ind(bf_j)) - c1*hxht*grad_v_t_eval(bf_i)*grad_tau_x_eval(bf_j)*qd_weights(pt_i)*qd_weights(pt_j);

                       % J_uu linear part: <v_t, v_t> + <v_x, v_x>
                       J_uu_lin(curr_ind(bf_i), curr_ind(bf_j)) = J_uu_lin(curr_ind(bf_i), curr_ind(bf_j)) + c1*hxht*grad_v_t_eval(bf_i)*grad_v_t_eval(bf_j)*qd_weights(pt_i)*qd_weights(pt_j);
                       J_uu_lin(curr_ind(bf_i), curr_ind(bf_j)) = J_uu_lin(curr_ind(bf_i), curr_ind(bf_j)) + c2*hxht*diff_const^2*grad_tau_x_eval(bf_i)*grad_tau_x_eval(bf_j)*qd_weights(pt_i)*qd_weights(pt_j);
                      
                       %%%%%%%%%% additional terms from a1*u %%%%%%%%%%%%%%
                       % J_su_init : < div(\tau), a1*v>
                       J_su_a1(curr_ind(bf_i), curr_ind(bf_j)) = J_su_a1(curr_ind(bf_i), curr_ind(bf_j)) + c1*hxht*grad_tau_x_eval(bf_i)*a1(qd_pt_x, qd_pt_t)*bf_eval(bf_j)*qd_weights(pt_i)*qd_weights(pt_j);
                       J_us_a1(curr_ind(bf_i), curr_ind(bf_j)) = J_us_a1(curr_ind(bf_i), curr_ind(bf_j)) + c1*hxht*grad_tau_x_eval(bf_j)*a1(qd_pt_x, qd_pt_t)*bf_eval(bf_i)*qd_weights(pt_i)*qd_weights(pt_j);

                       % J_uu_init : - <u_t, a1*v> - <v_t, a1*u > + <a1*u, a1*v>
                       J_uu_a1(curr_ind(bf_i), curr_ind(bf_j)) = J_uu_a1(curr_ind(bf_i), curr_ind(bf_j)) - c1*hxht*grad_v_t_eval(bf_j)*a1(qd_pt_x, qd_pt_t)*bf_eval(bf_i)*qd_weights(pt_i)*qd_weights(pt_j);
                       J_uu_a1(curr_ind(bf_i), curr_ind(bf_j)) = J_uu_a1(curr_ind(bf_i), curr_ind(bf_j)) - c1*hxht*grad_v_t_eval(bf_i)*a1(qd_pt_x, qd_pt_t)*bf_eval(bf_j)*qd_weights(pt_i)*qd_weights(pt_j);
                       J_uu_a1(curr_ind(bf_i), curr_ind(bf_j)) = J_uu_a1(curr_ind(bf_i), curr_ind(bf_j)) + c1*hxht*a1(qd_pt_x, qd_pt_t)*bf_eval(bf_j)*a1(qd_pt_x, qd_pt_t)*bf_eval(bf_i)*qd_weights(pt_i)*qd_weights(pt_j);

                    end
                end
                   
            end
        end
     end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% in order to improve initial guess, solve for linear part of f

J_ss_init = J_ss_lin;
J_su_init = J_su_lin + J_su_a1;
J_us_init = J_us_lin + J_us_a1;

J_uu_init = J_uu_lin + J_uu_a1;

%%%% impose boundary conditions on J_uu_init

J_uu_init(1:Nx_pts,:) = 0;
J_uu_init(1:Nx_pts:end) = 0;
J_uu_init(Nx_pts:Nx_pts:end) = 0;

J_us_init(1:Nx_pts,:) = 0;
J_us_init(1:Nx_pts:end) = 0;
J_us_init(Nx_pts:Nx_pts:end) = 0;

J_uu_init(1:Nx_pts, 1:Nx_pts) = eye(Nx_pts, Nx_pts);
J_uu_init(1:Nx_pts:end, 1:Nx_pts:end) = eye(Nt_pts, Nt_pts);
J_uu_init(Nx_pts:Nx_pts:end, Nx_pts:Nx_pts:end) = eye(Nt_pts, Nt_pts);

% ASSEMBLE MATRIX 
J_init_guess = [J_ss_init, J_su_init; J_us_init, J_uu_init];
block_size = size(J_ss_init,1);

% RHS 
F_init = [J_s_a0; J_u_a0];
F_init(block_size+1:block_size+Nx_pts) = u0;
F_init(block_size+1:Nx_pts:end) = bdy_left;
F_init(block_size+Nx_pts:Nx_pts:end) = bdy_right;

%%%%%%%%%%%%%%%%%%%%% SOLVE  %%%%%%%%%%%%%%%%%%%%%%
sol = J_init_guess \ F_init;

% INITIAL GUESS
sigma = sol(1:block_size);
u = sol(block_size+1:end);

% reshape
u_mat_init = zeros(Nx_pts, Nt_pts);

for ts=1:Nt_pts
    u_mat_init(:, ts) = u((ts-1)*Nx_pts+1:ts*Nx_pts);
end

figure
mesh(x, t(1:end), transpose(u_mat_init(:,1:end)));
xlabel('space');
ylabel('time');
zlabel('u');
%zlim([min(u), max(u)]);
title('heat equation LSFEM nonlinear');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% NEWTON ITERATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% reconstruckting NONLINEAR PART in every iteration
grad_J = 10;

for iter=1:max_iter
    if(and(norm(grad_J) < eps, iter > 1))
        % fprintf('norm gradient J: %d, iteration: %d\n ', norm(grad_J), iter);
        break;
    end
    
% ARRAY ALLOCATION 

% belonging to gradient
tau_x_f_int_vec = zeros(tot_pts, 1);
sigma_x_df_int_vec = zeros(tot_pts,1);
f_df_int_vec = zeros(tot_pts,1);
v_t_f_int_vec = zeros(tot_pts,1);
u_t_df_int_vec = zeros(tot_pts,1);

% for non linear parts of J_us, J_suX
tau_x_df_int = zeros(tot_pts);

% for non linear parts of J_uu
df_df_int = zeros(tot_pts);
f_d2f_int = zeros(tot_pts);
u_t_d2f_int = zeros(tot_pts);
s_x_d2f_int = zeros(tot_pts);
v_t_df_int_1 = zeros(tot_pts);
v_t_df_int_2 = zeros(tot_pts);


%%%% START INNER LOOPS %%%%
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
                
                d2f_local = d2f_dep_xt(coeff_vec_u,qd_pt_x,qd_pt_t); 
                df_local = df_dep_xt(coeff_vec_u,qd_pt_x,qd_pt_t);
                f_local = f_dep_xt(coeff_vec_u,qd_pt_x,qd_pt_t);
                
                sigma_x_local = dg_x(coeff_vec_sigma, qd_pt_t);        % to compute sigma_x term
                u_t_local = dg_t(coeff_vec_u, qd_pt_x);
                
                % iterate though basis functions for each element

                for bf_j=1:4
                    
                    % NON LINEAR BITS FOR GRADIENT
                    % compute: <\tau_x, f(u)> 
                    tau_x_f_int_vec(curr_ind(bf_j)) = tau_x_f_int_vec(curr_ind(bf_j)) ...
                                           + c1*hxht*grad_tau_x_eval(bf_j)*f_local*qd_weights(pt_i)*qd_weights(pt_j);
                       
                    % compute: < f(u), f'(u) v > 
                    f_df_int_vec(curr_ind(bf_j)) = f_df_int_vec(curr_ind(bf_j)) ...
                                           + c1*hxht*bf_eval(bf_j)*f_local*df_local*qd_weights(pt_i)*qd_weights(pt_j);
                       
                    % compute: < \sigma_x, f'(u) v >
                    sigma_x_df_int_vec(curr_ind(bf_j)) = sigma_x_df_int_vec(curr_ind(bf_j)) ...
                                           + c1*hxht*bf_eval(bf_j)*df_local*sigma_x_local*qd_weights(pt_i)*qd_weights(pt_j);
                    
                    % compute: - < u_t, f'(u) v >
                    u_t_df_int_vec(curr_ind(bf_j)) = u_t_df_int_vec(curr_ind(bf_j)) ...
                                           - c1*hxht*bf_eval(bf_j)*df_local*u_t_local*qd_weights(pt_i)*qd_weights(pt_j);
                    
                    % compute: - < v_t, f(u) > 
                    v_t_f_int_vec(curr_ind(bf_j)) = v_t_f_int_vec(curr_ind(bf_j)) ...
                                           - c1*hxht*grad_v_t_eval(bf_j)*f_local*qd_weights(pt_i)*qd_weights(pt_j);
                   
                    for bf_i=1:4                                      
                     
                       % NON LINEAR BITS FOR HESSIAN
                       % compute <\tau_x, f'(u)v >
                       tau_x_df_int(curr_ind(bf_i), curr_ind(bf_j)) = tau_x_df_int(curr_ind(bf_i), curr_ind(bf_j)) ...
                                           + c1*hxht*grad_tau_x_eval(bf_i)*df_local*bf_eval(bf_j)*qd_weights(pt_i)*qd_weights(pt_j);
                       
                       % compute: <f'(u)w , f'(u)v>
                       df_df_int(curr_ind(bf_i), curr_ind(bf_j)) = df_df_int(curr_ind(bf_i), curr_ind(bf_j)) ...
                                          + c1*hxht*bf_eval(bf_i)*bf_eval(bf_j)*df_local*df_local*qd_weights(pt_i)*qd_weights(pt_j);
                       
                       % compute <f(u), w^t f''(u) v >
                       f_d2f_int(curr_ind(bf_i), curr_ind(bf_j)) = f_d2f_int(curr_ind(bf_i), curr_ind(bf_j)) ...
                                          + c1*hxht*bf_eval(bf_i)*bf_eval(bf_j)*f_local*d2f_local*qd_weights(pt_i)*qd_weights(pt_j);
                       
                       % compute - <u_t, w^t f''(u) v >
                       u_t_d2f_int(curr_ind(bf_i), curr_ind(bf_j)) = u_t_d2f_int(curr_ind(bf_i), curr_ind(bf_j)) ...
                                          - c1*hxht*bf_eval(bf_i)*bf_eval(bf_j)*u_t_local*d2f_local*qd_weights(pt_i)*qd_weights(pt_j);
                        
                       % compute <sigma_x, w^t f''(u) v>
                       s_x_d2f_int(curr_ind(bf_i), curr_ind(bf_j)) = s_x_d2f_int(curr_ind(bf_i), curr_ind(bf_j)) ...
                                          + c1*hxht*bf_eval(bf_i)*bf_eval(bf_j)*sigma_x_local*d2f_local*qd_weights(pt_i)*qd_weights(pt_j); 

                       % compute: - < w_t, f'(u) v > 
                       v_t_df_int_1(curr_ind(bf_i), curr_ind(bf_j)) = v_t_df_int_1(curr_ind(bf_i), curr_ind(bf_j)) ...
                                          - c1*hxht*grad_v_t_eval(bf_i)*df_local*bf_eval(bf_j)*qd_weights(pt_i)*qd_weights(pt_j);
                       
                       % compute: - < v_t, f'(u) w >              
                       v_t_df_int_2(curr_ind(bf_i), curr_ind(bf_j)) = v_t_df_int_2(curr_ind(bf_i), curr_ind(bf_j)) ...
                                          - c1*hxht*grad_v_t_eval(bf_j)*df_local*bf_eval(bf_i)*qd_weights(pt_i)*qd_weights(pt_j);
                    end
                end
                   
            end
        end
     end
    
end

%%%%%%%%%%%%%%%%%% putting things together: gradient J

J_s = J_ss_lin*sigma + J_su_lin*u + tau_x_f_int_vec;
J_u = J_us_lin*sigma + J_uu_lin*u + f_df_int_vec + sigma_x_df_int_vec + v_t_f_int_vec + u_t_df_int_vec;

J_u(1:Nx_pts) = 0;
J_u(1:Nx_pts:end) = 0;
J_u(Nx_pts:Nx_pts:end) = 0;

grad_J = [J_s; J_u];
%norm(J_s)
%norm(J_u)
fprintf('norm gradient J: %d, iteration: %d\n ', norm(grad_J), iter);
%%%%%%%%%%%%%%%%%% putting things together: HESSIAN
J_ss = J_ss_lin;
J_uu = J_uu_lin + df_df_int + f_d2f_int + u_t_d2f_int + s_x_d2f_int + v_t_df_int_1 + v_t_df_int_2;

J_us = J_us_lin + tau_x_df_int';
J_su = J_su_lin + tau_x_df_int;

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

% assembling the blocks
Hess_J = [J_ss, J_su; J_us, J_uu];

%%%%%%%%%%%%%%%%%%%%%%% SOLVE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sol = sol - omega*Hess_J \ grad_J;

sigma = sol(1:block_size);
u = sol(block_size+1:end);

% enforcing bdy conditions on u
u(1:Nx_pts) = u0;
u(1:Nx_pts:end) = bdy_left;
u(Nx_pts:Nx_pts:end) = bdy_right;

sol = [sigma; u];

end
%% PLOT 

% exact solution
% figure
% mesh(x, t(1:end), transpose(u_exact(:,1:end)));
% xlabel('space');
% ylabel('time');
% zlabel('u');
% %zlim([min(u_exact), max(u_exact)]);
% title('exact solution');

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
%zlim([min(u), max(u)]);
title('heat equation LSFEM nonlinear');

norm(u_exact - u_mat)
norm(u_exact - u_mat_init)