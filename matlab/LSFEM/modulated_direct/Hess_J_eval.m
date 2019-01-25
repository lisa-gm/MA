% ******************************************************************* %
% ******************** compute Hessian J **************************** %
% ******************************************************************* %

function Hess_J = Hess_J_eval(hx, Nx_elem, ht, Nt_elem, bdy_cond, sigma, u, J_ss_lin, J_su_lin, J_us_lin, J_uu_lin, c1)
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


% ARRAY ALLOCATION 
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
    
     coeff_vec_sigma = sigma(curr_ind);
     coeff_vec_u = u(curr_ind);
     
     [coeff_vec_f,coeff_vec_df,coeff_vec_d2f] = F_eval(coeff_vec_u);
     
        for pt_j=1:qd_deg

            for pt_i=1:qd_deg             
                
                bf_eval = basis_fct_eval(qd_pts(pt_i),qd_pts(pt_j));   % eval basis fct in quad pts       
                grad_tau_x_eval = grad_x_fct_eval(qd_pts(pt_i),qd_pts(pt_j));
                grad_v_t_eval = grad_t_fct_eval(qd_pts(pt_i),qd_pts(pt_j));              
                
                % iterate though basis functions for each element
                
                for bf_j=1:4
                    for bf_i=1:4                        
                       
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

end