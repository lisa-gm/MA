% ******************************************************** %
% ***************** NEWTON STEP ************************** %
% ******************************************************** %


function [sigma, u, norm_grad_J] = newton_step_lsfem(S, Nx_elem, T, Nt_elem, u0, bdy_cond, bdy_left, bdy_right, sigma, u, a0, a1, a2, a3, alpha, J_ss_lin, J_su_lin, J_us_lin, J_uu_lin, c1, c2, diff_const, eps, iter)
max_iter = 50;
damp_fact = 1/3;
sol = [sigma; u];

int_space = [0,S];
hx = (int_space(2) - int_space(1))/Nx_elem;
Nx_pts = Nx_elem+1;

int_time = [0,T];
ht = (int_time(2) - int_time(1))/Nt_elem;
Nt_pts = Nt_elem+1;

tot_pts = Nx_pts*Nt_pts;

hxht = ht*hx;

%% compute indices of boundary IDs

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
% f = a0 + a1*u + a2*u^2 + a3*u^3
% a0 = @(x,t) 0; a1 = @(x,t) 0; a2 = @(x,t) 0; a3 = 0;

% approximate for integration, should be exact for polynomials of deg 1
d2f_dep_xt = @(coeff_vec,x,t, pt_i, pt_j) alpha*(2*a2(x,t) + 6*a3*dot(coeff_vec, basis_fct_eval(qd_pts(pt_i), qd_pts(pt_j))));
df_dep_xt = @(coeff_vec,x,t, pt_i, pt_j) alpha*(a1(x,t) + 2*a2(x,t)*dot(coeff_vec, basis_fct_eval(qd_pts(pt_i), qd_pts(pt_j))) + 3*a3*(dot(coeff_vec, basis_fct_eval(qd_pts(pt_i), qd_pts(pt_j))))^2);
f_dep_xt = @(coeff_vec,x,t, pt_i, pt_j) alpha*(a0(x,t) + a1(x,t)*dot(coeff_vec, basis_fct_eval(qd_pts(pt_i), qd_pts(pt_j))) + ...
            a2(x,t)*(dot(coeff_vec, basis_fct_eval(qd_pts(pt_i), qd_pts(pt_j))))^2 + a3*(dot(coeff_vec, basis_fct_eval(qd_pts(pt_i), qd_pts(pt_j))))^3);

% ARRAY ALLOCATION 

% belonging to gradient
tau_x_f_int_vec = zeros(tot_pts, 1);
sigma_x_df_int_vec = zeros(tot_pts,1);
f_df_int_vec = zeros(tot_pts,1);
v_t_f_int_vec = zeros(tot_pts,1);
u_t_df_int_vec = zeros(tot_pts,1);

% for non linear parts of J_us, J_su
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
                     
                       % external residual computation
                       %res_J_s(curr_ind(bf_j)) = res_J_s(curr_ind(bf_j)) + c2*hxht*bf_eval(bf_j)*coeff_vec_sigma(bf_j)*bf_eval(bf_i)*qd_weights(pt_i)*qd_weights(pt_j); 
                        
                        
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

%norm(J_u - J_u_test)
%norm(J_s - J_s_test)

J_u(1:Nx_pts) = 0;

if(strcmp(bdy_cond, 'Dirichlet'))
    J_u(1:Nx_pts:end) = 0;
    J_u(Nx_pts:Nx_pts:end) = 0;
end

if(strcmp(bdy_cond, 'Neumann'))
    J_s(1:Nx_pts:end) = 0;
    J_s(Nx_pts:Nx_pts:end) = 0;
end

grad_J = [J_s; J_u];

%norm(J_s)
%norm(J_u)
norm_grad_J = norm(grad_J);

fprintf('norm residual: %d, iteration: %d\n ', norm_grad_J, iter);

if(norm_grad_J > eps || iter == 1)
    
%%%%%%%%%%%%%%%%%% putting things together: HESSIAN
J_ss = J_ss_lin;
J_uu = J_uu_lin + df_df_int + f_d2f_int + u_t_d2f_int + s_x_d2f_int + v_t_df_int_1 + v_t_df_int_2;

J_us = J_us_lin + tau_x_df_int';
J_su = J_su_lin + tau_x_df_int;

%%%%%%%%%%%% putting in Dirichlet boundary conditions

J_uu(1:Nx_pts,:) = 0;
J_us(1:Nx_pts,:) = 0;
J_uu(1:Nx_pts, 1:Nx_pts) = eye(Nx_pts, Nx_pts);


if(strcmp(bdy_cond, 'Dirichlet'))
    J_uu(1:Nx_pts:end) = 0;
    J_uu(Nx_pts:Nx_pts:end) = 0;

    J_us(1:Nx_pts:end) = 0;
    J_us(Nx_pts:Nx_pts:end) = 0;

    J_uu(1:Nx_pts:end, 1:Nx_pts:end) = eye(Nt_pts, Nt_pts);
    J_uu(Nx_pts:Nx_pts:end, Nx_pts:Nx_pts:end) = eye(Nt_pts, Nt_pts);
end


if(strcmp(bdy_cond, 'Neumann'))
    J_ss(1:Nx_pts:end) = 0;
    J_ss(Nx_pts:Nx_pts:end) = 0;

    J_su(1:Nx_pts:end) = 0;
    J_su(Nx_pts:Nx_pts:end) = 0;

    J_ss(1:Nx_pts:end, 1:Nx_pts:end) = eye(Nt_pts, Nt_pts);
    J_ss(Nx_pts:Nx_pts:end, Nx_pts:Nx_pts:end) = eye(Nt_pts, Nt_pts);
end

% assembling the blocks
Hess_J = [J_ss, J_su; J_us, J_uu];

my_u_old = u;
sigma_old = sigma;
%%%%%%%%%%%%%%%%%%%%%%% SOLVE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
update = - Hess_J \ grad_J;

%%%%%%%%%%%%%%%%%%%%%% CHECK UPDATE %%%%%%%%%%%%%%%%%%%%%%%%%%%
temp = dot(update([inner_ind_sigma, tot_pts+inner_ind_u]), grad_J([inner_ind_sigma, tot_pts+inner_ind_u]));
fprintf(' sign(dot( update, grad J)) : %d \n', sign(temp));
if(temp > 0)
    % use opposite direction since problem not convex in current direction
    update = - update;
    fprintf('non convex \n');
end

% compute J of current solution to compare with potential updates
J_eval_old = J_eval_lsfem(S, Nx_elem, T, Nt_elem, sigma, u, a0, a1, a2, a3, alpha, bdy_cond, c1, c2, diff_const);

%sol = sol + update;
temp_damp_fact = damp_fact;
inner_it = 0;

% sigma_temp = sigma;
% u_temp = u;

temp_sol = sol;

while(inner_it < 50)    
    
     %sigma_temp(inner_ind_sigma) = u(inner_ind_sigma) + temp_damp_fact*update(inner_ind_sigma);
     %u_temp(inner_ind_u) = u(inner_ind_u) + temp_damp_fact*update(tot_pts + inner_ind_u);
    temp_sol([inner_ind_sigma, tot_pts+inner_ind_u]) = temp_sol([inner_ind_sigma, tot_pts+inner_ind_u]) ...
                        + damp_fact*update([inner_ind_sigma, tot_pts+inner_ind_u]);

    sigma_temp = temp_sol(1:tot_pts);    
    u_temp = temp_sol(tot_pts+1:end);
    
    J_eval_temp = J_eval_lsfem(S, Nx_elem, T, Nt_elem, sigma_temp, u_temp, a0, a1, a2, a3, alpha, bdy_cond, c1, c2, diff_const);

    fprintf('damping factor = %g, new J = %g, old J = %g\n',temp_damp_fact,J_eval_temp,J_eval_old);
    
    if J_eval_temp < J_eval_old
      break;
    end
    
    temp_damp_fact = temp_damp_fact/2;
    inner_it = inner_it + 1;
    if(inner_it == 50)
       fprintf('max iter reached!');
    end
end
  
  sigma = sigma_temp;
  u = u_temp;
  
end
   
% enforcing bdy conditions on u
u(1:Nx_pts) = u0;

    if(strcmp(bdy_cond, 'Dirichlet'))
        u(1:Nx_pts:end) = bdy_left;
        u(Nx_pts:Nx_pts:end) = bdy_right;
    end

    if(strcmp(bdy_cond, 'Neumann'))
        sigma(1:Nx_pts:end) = -bdy_left;
        sigma(Nx_pts:Nx_pts:end) = bdy_right;
    end
my_u = u;
Hess_J_no_bdy = Hess_J([inner_ind_sigma, tot_pts+inner_ind_u], [inner_ind_sigma, tot_pts+inner_ind_u]);

J_ss_no_bdy = J_ss(inner_ind_sigma, inner_ind_sigma);
J_su_no_bdy = J_su(inner_ind_sigma, inner_ind_u);
J_us_no_bdy = J_us(inner_ind_u, inner_ind_sigma);
J_uu_no_bdy = J_uu(inner_ind_u, inner_ind_u);


J_s_no_bdy = grad_J(inner_ind_sigma);
J_u_no_bdy = grad_J(inner_ind_u);
save('my_variables.mat', 'update', 'J_ss_no_bdy', 'J_su_no_bdy', 'J_us_no_bdy', 'J_uu_no_bdy', 'Hess_J_no_bdy', 'J_s_no_bdy', 'J_u_no_bdy', 'J_eval_temp', 'sigma', 'my_u', 'sigma_old', 'my_u_old');

sol = [sigma; u];
end