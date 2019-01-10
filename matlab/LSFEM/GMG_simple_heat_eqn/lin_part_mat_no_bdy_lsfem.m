% ******************************************************** %
% *************** SET UP LINEAR MATRICES ***************** %
% ******************************************************** %

function [J_ss_lin, J_su_lin, J_us_lin, J_uu_lin] = lin_part_mat_no_bdy_lsfem (S, Nx_elem, T, Nt_elem, c1 ,c2, diff_const)

int_space = [0,S];
hx = (int_space(2) - int_space(1))/Nx_elem;
Nx_pts = Nx_elem+1;

int_time = [0,T];
ht = (int_time(2) - int_time(1))/Nt_elem;
Nt_pts = Nt_elem+1;

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
grad_x_fct_eval = @(x,t) S/hx*[(-1)*(1-t); 1*(1-t); (-1)*t;  1*t];
grad_t_fct_eval = @(x,t) T/ht*[(1-x)*(-1); x*(-1); (1-x)*1;  x*1];


% linear parts of Hessian
J_ss_lin = zeros(tot_pts);
J_su_lin = zeros(tot_pts);                          % compute both to check symmetry
J_us_lin = zeros(tot_pts);
J_uu_lin = zeros(tot_pts);
        
%% loop construct matrices

for elem_j = 1:Nt_elem                                      % iterate through all elements
    for elem_i = 1:Nx_elem
        
     curr_ind = get_elem_ind(elem_i, elem_j);               % get vector with indices of current elem
          
        for pt_j=1:qd_deg                                   % iterate through quadrature points 

            for pt_i=1:qd_deg             
                
                bf_eval = basis_fct_eval(qd_pts(pt_i),qd_pts(pt_j));   % eval basis fct in quad pts       
                grad_tau_x_eval = grad_x_fct_eval(qd_pts(pt_i),qd_pts(pt_j));
                grad_v_t_eval = grad_t_fct_eval(qd_pts(pt_i),qd_pts(pt_j));
                                                       
                for bf_j=1:4                                           % iterate though basis functions
                              
                    for bf_i=1:4
                       % J_ss linear part: <\tau, \tau> + <\tau_x, \tau_x> 
                       J_ss_lin(curr_ind(bf_i), curr_ind(bf_j)) = J_ss_lin(curr_ind(bf_i), curr_ind(bf_j)) + c2*hxht*bf_eval(bf_i)*bf_eval(bf_j)*qd_weights(pt_i)*qd_weights(pt_j);
                       J_ss_lin(curr_ind(bf_i), curr_ind(bf_j)) = J_ss_lin(curr_ind(bf_i), curr_ind(bf_j)) + c1*hxht*grad_tau_x_eval(bf_i)*grad_tau_x_eval(bf_j)*qd_weights(pt_i)*qd_weights(pt_j);
                       
                       % J_su linear part: -diff_const < \tau, v_x> - <v_t, \tau_x>
                       J_su_lin(curr_ind(bf_i), curr_ind(bf_j)) = J_su_lin(curr_ind(bf_i), curr_ind(bf_j)) - c2*hxht*diff_const*bf_eval(bf_i)*grad_tau_x_eval(bf_j)*qd_weights(pt_i)*qd_weights(pt_j);
                       J_su_lin(curr_ind(bf_i), curr_ind(bf_j)) = J_su_lin(curr_ind(bf_i), curr_ind(bf_j)) - c1*hxht*grad_v_t_eval(bf_j)*grad_tau_x_eval(bf_i)*qd_weights(pt_i)*qd_weights(pt_j);
                       
                       % J_us linear part
                       J_us_lin(curr_ind(bf_i), curr_ind(bf_j)) = J_us_lin(curr_ind(bf_i), curr_ind(bf_j)) - c2*hxht*diff_const*bf_eval(bf_j)*grad_tau_x_eval(bf_i)*qd_weights(pt_i)*qd_weights(pt_j);
                       J_us_lin(curr_ind(bf_i), curr_ind(bf_j)) = J_us_lin(curr_ind(bf_i), curr_ind(bf_j)) - c1*hxht*grad_v_t_eval(bf_i)*grad_tau_x_eval(bf_j)*qd_weights(pt_i)*qd_weights(pt_j);

                       % J_uu linear part: <v_t, v_t> + <v_x, v_x>
                       J_uu_lin(curr_ind(bf_i), curr_ind(bf_j)) = J_uu_lin(curr_ind(bf_i), curr_ind(bf_j)) + c1*hxht*grad_v_t_eval(bf_i)*grad_v_t_eval(bf_j)*qd_weights(pt_i)*qd_weights(pt_j);
                       J_uu_lin(curr_ind(bf_i), curr_ind(bf_j)) = J_uu_lin(curr_ind(bf_i), curr_ind(bf_j)) + c2*hxht*diff_const^2*grad_tau_x_eval(bf_i)*grad_tau_x_eval(bf_j)*qd_weights(pt_i)*qd_weights(pt_j);
                      
                     
                    end
                end
                   
            end
        end
    end   
end

end