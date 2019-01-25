% ***************************************************************** %
% ****************** compute gradient J at u ********************** %
% ***************************************************************** %

function grad_J = grad_J_eval(hx, Nx_elem, ht, Nt_elem, bdy_cond, sigma, u, J_ss_lin, J_su_lin, J_us_lin, J_uu_lin, c1)
%%%%%%%%%%%%%%% SET UP %%%%%%%%%%%%%%
Nx_pts = Nx_elem + 1;
Nt_pts = Nt_elem + 1;
tot_pts = Nx_pts*Nt_pts;

hxht = ht*hx;

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

%%%% START INNER LOOPS %%%%
for elem_j = 1:Nt_elem                                      % iterate through all elements
    for elem_i = 1:Nx_elem
        
     curr_ind = get_elem_ind(elem_i, elem_j);               % get vector with indices of current elem
     
     coeff_vec_sigma = sigma(curr_ind);
     coeff_vec_u = u(curr_ind);
     
     [coeff_vec_f,coeff_vec_df,~] = F_eval(coeff_vec_u);
     
        for pt_j=1:qd_deg

            for pt_i=1:qd_deg             
                
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
                       
                    end
                end
                   
            end
            
        end
     end
    
end

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

end