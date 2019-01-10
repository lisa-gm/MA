% ********************************************************************* %
% ************* CONSTRUCT RHS OF POISSON PROBLEM ********************** %
% ********************************************************************* %

% rhs = - <v_t, f> + < div(\sigma), f>

function rhs = poisson_rhs_no_bdy_lsfem (S, Nx_elem, f)
int_space = [0,S];
hx = (int_space(2) - int_space(1))/Nx_elem;
Nx_pts = Nx_elem+1;

f_u = f(Nx_pts+1: end);

%% quadrature
qd_deg = 3;
qd_pts = 0.5*([-sqrt(3/5),0,+sqrt(3/5)]+1);                         %3-point gauss on [0,1]
qd_weights = 0.5*[5/9,8/9,5/9];

%% basis functions, and their gradients

basis_fct_eval = @(x) [(1-x); x];
grad_x_fct_eval = @(x) 1/hx*[(-1); 1];

% linear parts of Hessian
rhs_sigma = zeros(Nx_pts, 1);
rhs_u = zeros(Nx_pts,1);
        
%% loop construct matrices

for elem_i = 1:Nx_elem
  
    curr_ind = [elem_i; elem_i+1];
    f_curr = f_u(curr_ind);
    
        for pt_i=1:qd_deg             
                
         bf_eval = basis_fct_eval(qd_pts(pt_i));   % eval basis fct in quad pts       
         grad_tau_x_eval = grad_x_fct_eval(qd_pts(pt_i));
                                                       
                for bf_j=1:2                                          % iterate though basis functions                             
                   for bf_i=1:2
                       
                       rhs_u(curr_ind(bf_i)) = rhs_u(curr_ind(bf_i)) + hx*grad_tau_x_eval(bf_i)*f_curr(bf_j)*bf_eval(bf_j)*qd_weights(pt_i);                   
                   
                   end
                end
                   
        end
end   

rhs = [rhs_sigma; rhs_u];

end