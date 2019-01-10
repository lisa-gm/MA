% ******************************************************** %
% *************** SET UP LINEAR MATRICES ***************** %
% ******************************************************** %

function [J_ss_lin, J_su_lin, J_us_lin, J_uu_lin] = poisson_mat_no_bdy_lsfem (S, Nx_elem, diff_const)

int_space = [0,S];
hx = (int_space(2) - int_space(1))/Nx_elem;
Nx_pts = Nx_elem+1;

%% quadrature
qd_deg = 3;
qd_pts = 0.5*([-sqrt(3/5),0,+sqrt(3/5)]+1);                         %3-point gauss on [0,1]
qd_weights = 0.5*[5/9,8/9,5/9];

%% basis functions, and their gradients

basis_fct_eval = @(x) [(1-x); x];
grad_x_fct_eval = @(x) S/hx*[(-1); 1];

% linear parts of Hessian
J_ss_lin = zeros(Nx_pts);
J_su_lin = zeros(Nx_pts);                          % compute both to check symmetry
J_us_lin = zeros(Nx_pts);
J_uu_lin = zeros(Nx_pts);
        
%% loop construct matrices

for elem_i = 1:Nx_elem
  
    curr_ind = [elem_i; elem_i+1];
        for pt_i=1:qd_deg             
                
         bf_eval = basis_fct_eval(qd_pts(pt_i));   % eval basis fct in quad pts       
         grad_tau_x_eval = grad_x_fct_eval(qd_pts(pt_i));
                                                       
                for bf_j=1:2                                          % iterate though basis functions                             
                   for bf_i=1:2
                       % J_ss linear part: <\tau, \tau> + <\tau_x, \tau_x> 
                       J_ss_lin(curr_ind(bf_i), curr_ind(bf_j)) = J_ss_lin(curr_ind(bf_i), curr_ind(bf_j)) + hx*bf_eval(bf_i)*bf_eval(bf_j)*qd_weights(pt_i);
                       J_ss_lin(curr_ind(bf_i), curr_ind(bf_j)) = J_ss_lin(curr_ind(bf_i), curr_ind(bf_j)) + hx*grad_tau_x_eval(bf_i)*grad_tau_x_eval(bf_j)*qd_weights(pt_i);
                       
                       % J_su linear part: - diff_const < \tau, v_x> 
                       J_su_lin(curr_ind(bf_i), curr_ind(bf_j)) = J_su_lin(curr_ind(bf_i), curr_ind(bf_j)) - hx*diff_const*bf_eval(bf_i)*grad_tau_x_eval(bf_j)*qd_weights(pt_i);
                       
                       % J_us linear part
                       J_us_lin(curr_ind(bf_i), curr_ind(bf_j)) = J_us_lin(curr_ind(bf_i), curr_ind(bf_j)) - hx*diff_const*bf_eval(bf_j)*grad_tau_x_eval(bf_i)*qd_weights(pt_i);

                       % J_uu linear part: <v_x, v_x>
                       J_uu_lin(curr_ind(bf_i), curr_ind(bf_j)) = J_uu_lin(curr_ind(bf_i), curr_ind(bf_j)) + hx*diff_const^2*grad_tau_x_eval(bf_i)*grad_tau_x_eval(bf_j)*qd_weights(pt_i);
                      
                     
                   end
                end
                   
        end
end   

end