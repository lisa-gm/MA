% ************************************************************ %
% ************* evaluate functional J ************************ %
% ************************************************************ %

function J_eval = J_eval_lsfem(S, Nx_elem, T, Nt_elem, sigma, u, c1, c2, diff_const)
% J= c1 || u_t - div(sigma) - f(u) ||^2+c2 ||sigma-diff_const*nabla u||^2

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
grad_x_fct_eval = @(x,t) 1/hx*[(-1)*(1-t); 1*(1-t); (-1)*t;  1*t];
grad_t_fct_eval = @(x,t) 1/ht*[(1-x)*(-1); x*(-1); (1-x)*1;  x*1];

% t1 = || u_t - div(sigma) - f(u) ||^2, t2 =||sigma-diff_const*nabla(u)||^2
% compute as sum

t1 = 0;
t2 = 0;

%%%% START INNER LOOPS %%%%
for elem_j = 1:Nt_elem                                      % iterate through all elements
    for elem_i = 1:Nx_elem
        
     curr_ind = get_elem_ind(elem_i, elem_j);               % get vector with indices of current elem
     
     curr_elem_x = int_space(1) + (elem_i-1)*hx;            % get global coordinate, lower left corner
     curr_elem_t = int_time(1) + (elem_j-1)*ht;             % of current elem
     
     coeff_vec_sigma = sigma(curr_ind);
     coeff_vec_u = u(curr_ind);
     
     %coeff_vec_u = [1;2;3;4];
     
     [coeff_vec_f,~, ~] = F_eval(coeff_vec_u);


     
        for pt_j=1:qd_deg

            for pt_i=1:qd_deg             
                
                bf_eval = basis_fct_eval(qd_pts(pt_i),qd_pts(pt_j));   % eval basis fct in quad pts        
                grad_tau_x_eval = grad_x_fct_eval(qd_pts(pt_i),qd_pts(pt_j));
                grad_v_t_eval = grad_t_fct_eval(qd_pts(pt_i),qd_pts(pt_j)); 
                
                % iterate though basis functions for each element

                for bf_j=1:4
                    for bf_i=1:4
                        
                t1 = t1 + 0.5*c1*hxht*(coeff_vec_u(bf_j)*grad_v_t_eval(bf_j) - coeff_vec_sigma(bf_j)*grad_tau_x_eval(bf_j) - coeff_vec_f(bf_j)*bf_eval(bf_j))...
                   *(coeff_vec_u(bf_i)*grad_v_t_eval(bf_i) - coeff_vec_sigma(bf_i)*grad_tau_x_eval(bf_i) - coeff_vec_f(bf_i)*bf_eval(bf_i))*qd_weights(pt_i)*qd_weights(pt_j);                    

                t2 = t2 + 0.5*c2*hxht*(coeff_vec_sigma(bf_j)*bf_eval(bf_j)-diff_const*coeff_vec_u(bf_j)*grad_tau_x_eval(bf_j))...
                       *(coeff_vec_sigma(bf_i)*bf_eval(bf_i)-diff_const*coeff_vec_u(bf_i)*grad_tau_x_eval(bf_i))*qd_weights(pt_i)*qd_weights(pt_j);
                                       
                    end
                end
                   
            end
        end
     end
    
end

J_eval = t1 + t2;


end

