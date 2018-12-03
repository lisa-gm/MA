% ******************************************************** %
% ***************** NEWTON STEP ************************** %
% ******************************************************** %


function [sigma, u, norm_grad_J, damp_fact] = newton_step_lsfem(S, Nx_elem, T, Nt_elem, u0, bdy_cond, bdy_left, bdy_right, sigma, u, J_ss_lin, J_su_lin, J_us_lin, J_uu_lin, c1, c2, diff_const, eps, damp_fact, iter)

%%%%%%%%%%%%%%% SET UP %%%%%%%%%%%%%%
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
%% basis functions, and their gradients

basis_fct_eval = @(x,t) [(1-x)*(1-t); x*(1-t); (1-x)*t;  x*t];
grad_x_fct_eval = @(x,t) 1/hx*[(-1)*(1-t); 1*(1-t); (-1)*t;  1*t];
grad_t_fct_eval = @(x,t) 1/ht*[(1-x)*(-1); x*(-1); (1-x)*1;  x*1];

% function F=basis_function(x,t)
% F = [(1-x)*(1-t); x*(1-t); (1-x)*t;  x*t];
% end

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
     
     [coeff_vec_f,coeff_vec_df,coeff_vec_d2f] = F_eval(coeff_vec_u);
     
        for pt_j=1:qd_deg
            qd_pt_t = curr_elem_t + qd_pts(pt_j)*ht;            % get global quad pt in y, scaled

            for pt_i=1:qd_deg             
                qd_pt_x = curr_elem_x + qd_pts(pt_i)*hx;        % get global quad pt in x, scaled
                
                bf_eval = basis_fct_eval(qd_pts(pt_i),qd_pts(pt_j));   % eval basis fct in quad pts       
                grad_tau_x_eval = grad_x_fct_eval(qd_pts(pt_i),qd_pts(pt_j));
                grad_v_t_eval = grad_t_fct_eval(qd_pts(pt_i),qd_pts(pt_j));              
                
                % iterate though basis functions for each element
                
                for bf_j=1:4
                    for bf_i=1:4 
                                                
                    % NON LINEAR BITS FOR GRADIENT
                    % compute: <\tau_x, f(u)> 
                    tau_x_f_int_vec(curr_ind(bf_j)) = tau_x_f_int_vec(curr_ind(bf_j)) ...
                                           + c1*hxht*grad_tau_x_eval(bf_j)*coeff_vec_f(bf_i)*bf_eval(bf_i)*qd_weights(pt_i)*qd_weights(pt_j);
                    
                    % compute: < f(u), f'(u) v > 
                    f_df_int_vec(curr_ind(bf_j)) = f_df_int_vec(curr_ind(bf_j)) ...
                                           + c1*hxht*coeff_vec_f(bf_i)*bf_eval(bf_i)*bf_eval(bf_j)*coeff_vec_df(bf_j)*qd_weights(pt_i)*qd_weights(pt_j);
                       
                    % compute: < \sigma_x, f'(u) v >
                    sigma_x_df_int_vec(curr_ind(bf_j)) = sigma_x_df_int_vec(curr_ind(bf_j)) ...
                                           + c1*hxht*coeff_vec_sigma(bf_i)*grad_tau_x_eval(bf_i)*bf_eval(bf_j)*coeff_vec_df(bf_j)*qd_weights(pt_i)*qd_weights(pt_j);
                    
                    % compute: - < u_t, f'(u) v >
                    u_t_df_int_vec(curr_ind(bf_j)) = u_t_df_int_vec(curr_ind(bf_j)) ...
                                           - c1*hxht*grad_v_t_eval(bf_i)*coeff_vec_u(bf_i)*bf_eval(bf_j)*coeff_vec_df(bf_j)*qd_weights(pt_i)*qd_weights(pt_j);
                    
                    % compute: - < v_t, f(u) > 
                    v_t_f_int_vec(curr_ind(bf_j)) = v_t_f_int_vec(curr_ind(bf_j)) ...
                                           - c1*hxht*grad_v_t_eval(bf_j)*coeff_vec_f(bf_i)*bf_eval(bf_i)*qd_weights(pt_i)*qd_weights(pt_j);                        
                       
                       
                       % NON LINEAR BITS FOR HESSIAN
                       % compute <\tau_x, f'(u)v >
                       tau_x_df_int(curr_ind(bf_i), curr_ind(bf_j)) = tau_x_df_int(curr_ind(bf_i), curr_ind(bf_j)) ...
                                           + c1*hxht*grad_tau_x_eval(bf_i)*coeff_vec_df(bf_j)*bf_eval(bf_j)*qd_weights(pt_i)*qd_weights(pt_j);
                       
                                       
                       % compute: <f'(u)w , f'(u)v>
                       df_df_int(curr_ind(bf_i), curr_ind(bf_j)) = df_df_int(curr_ind(bf_i), curr_ind(bf_j)) ...
                                          + c1*hxht*bf_eval(bf_i)*bf_eval(bf_j)*coeff_vec_df(bf_j)*coeff_vec_df(bf_i)*qd_weights(pt_i)*qd_weights(pt_j);
                       
                              
                       % compute <f(u), w^t f''(u) v >
                       f_d2f_int(curr_ind(bf_i), curr_ind(bf_i)) = f_d2f_int(curr_ind(bf_i), curr_ind(bf_i)) ...
                                          + c1*hxht*bf_eval(bf_i)*coeff_vec_d2f(bf_i)*bf_eval(bf_j)*coeff_vec_f(bf_j)*qd_weights(pt_i)*qd_weights(pt_j);
                       
                       % compute - <u_t, w^t f''(u) v >
                       u_t_d2f_int(curr_ind(bf_i), curr_ind(bf_i)) = u_t_d2f_int(curr_ind(bf_i), curr_ind(bf_i)) ...
                                           - c1*hxht*bf_eval(bf_i)*coeff_vec_d2f(bf_i)*coeff_vec_u(bf_j)*grad_v_t_eval(bf_j)*qd_weights(pt_i)*qd_weights(pt_j);
                       
                       % compute <sigma_x, w^t f''(u) v>
                       s_x_d2f_int(curr_ind(bf_i), curr_ind(bf_i)) = s_x_d2f_int(curr_ind(bf_i), curr_ind(bf_i)) ...
                                          + c1*hxht*bf_eval(bf_i)*coeff_vec_d2f(bf_i)*coeff_vec_sigma(bf_j)*grad_tau_x_eval(bf_j)*qd_weights(pt_i)*qd_weights(pt_j);
                                 
                       % compute: - < w_t, f'(u) v > 
                       v_t_df_int_1(curr_ind(bf_i), curr_ind(bf_j)) = v_t_df_int_1(curr_ind(bf_i), curr_ind(bf_j)) ...
                                          - c1*hxht*grad_v_t_eval(bf_i)*coeff_vec_df(bf_j)*bf_eval(bf_j)*qd_weights(pt_i)*qd_weights(pt_j);
                       
                       % compute: - < v_t, f'(u) w >              
                       v_t_df_int_2(curr_ind(bf_i), curr_ind(bf_j)) = v_t_df_int_2(curr_ind(bf_i), curr_ind(bf_j)) ...
                                          - c1*hxht*grad_v_t_eval(bf_j)*coeff_vec_df(bf_i)*bf_eval(bf_i)*qd_weights(pt_i)*qd_weights(pt_j);
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

if(strcmp(bdy_cond, 'Dirichlet'))
    J_u(1:Nx_pts:end) = 0;
    J_u(Nx_pts:Nx_pts:end) = 0;
end

if(strcmp(bdy_cond, 'Neumann'))
    J_s(1:Nx_pts:end) = 0;
    J_s(Nx_pts:Nx_pts:end) = 0;
end

grad_J = [J_s; J_u];
norm_grad_J = norm(grad_J);

%%%%%%%%%%%%%%%%%% putting things together: HESSIAN %%%%%%%%%%%%%%%%%%%
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

Hess_J_inner = Hess_J([inner_ind_sigma, tot_pts+inner_ind_u], [inner_ind_sigma, tot_pts+inner_ind_u]);
grad_J_inner = grad_J([inner_ind_sigma, tot_pts+inner_ind_u]);

%%%%%%%%%%%%%%%%%%%%%%% SOLVE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
update_inner = - Hess_J_inner \ grad_J_inner;

%%%%%%%%%%%%%%%%%%%%%% CHECK UPDATE %%%%%%%%%%%%%%%%%%%%%%%%%%%
temp = dot(update_inner, grad_J_inner);
%fprintf(' sign(dot( update, grad J)) : %d \n', sign(temp));
if(temp > 0)
    % use opposite direction since problem not convex in current direction
    update_inner = - update_inner;
    fprintf('not convex \n');
end

% compute J of current solution to compare with potential updates
J_eval_old = J_eval_lsfem(S, Nx_elem, T, Nt_elem, sigma, u, c1, c2, diff_const);

inner_it = 0;

u_temp = u;
sigma_temp = sigma;

damp_fact = 1;
while(inner_it < 60)    
    
    sigma_temp(inner_ind_sigma) = sigma(inner_ind_sigma)...
            + damp_fact*update_inner(1:length(inner_ind_sigma));
    u_temp(inner_ind_u) = u(inner_ind_u) ...
            + damp_fact*update_inner(length(inner_ind_sigma)+1:end);
    
    J_eval_temp = J_eval_lsfem(S, Nx_elem, T, Nt_elem, sigma_temp, u_temp, c1, c2, diff_const);

    fprintf('damping factor = %g, new J = %g, old J = %g\n', damp_fact,J_eval_temp,J_eval_old);
    
    if J_eval_temp < J_eval_old
      break;
    end
    
    damp_fact = damp_fact/2;
    inner_it = inner_it + 1;
    if(inner_it == 50)
       fprintf('max iter reached!');
    end
end
  
  sigma = sigma_temp;
  u = u_temp;
     
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
   
end
