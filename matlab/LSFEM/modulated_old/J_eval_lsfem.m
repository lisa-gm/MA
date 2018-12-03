% ************************************************************ %
% ************* evaluate functional J ************************ %
% ************************************************************ %

function J_eval = J_eval_lsfem(S, Nx_elem, T, Nt_elem, sigma, u, a0, a1, a2, a3, alpha, bdy_cond, c1, c2, diff_const)
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

%% directional derivatives in direction x and t

dg_x = @(coeff_vec, t) 1/hx*(coeff_vec(1)*(-1)*(1-t) + coeff_vec(2)*(1-t) + coeff_vec(3)*t*(-1) + coeff_vec(4)*t);
dg_t = @(coeff_vec, x) 1/ht*(coeff_vec(1)*(1-x)*(-1) + coeff_vec(2)*x*(-1) + coeff_vec(3)*(1-x) + coeff_vec(4)*x);

%% basis functions, and their gradients

basis_fct_eval = @(x,t) [(1-x)*(1-t); x*(1-t); (1-x)*t;  x*t];

% function F=basis_function(x,t)
% F = [(1-x)*(1-t); x*(1-t); (1-x)*t;  x*t];
% end

%% define f, polynomial of the form: sum_i a_i u^i with a_i = a(x,t), i<=3
% f = a0 + a1*u + a2*u^2 + a3*u^3
% a0 = @(x,t) 0; a1 = @(x,t) 0; a2 = @(x,t) 0; a3 = 0;

% approximate for integration, should be exact for polynomials of deg 1
f_dep_xt = @(coeff_vec,x,t, pt_i, pt_j) alpha*(a0(x,t) + a1(x,t)*dot(coeff_vec, basis_fct_eval(qd_pts(pt_i), qd_pts(pt_j))) + ...
            a2(x,t)*(dot(coeff_vec, basis_fct_eval(qd_pts(pt_i), qd_pts(pt_j))))^2 + a3*(dot(coeff_vec, basis_fct_eval(qd_pts(pt_i), qd_pts(pt_j))))^3);

% t1 = || u_t - div(sigma) - f(u) ||^2
% compute as vector, sum afterwards

t1_int_vec = zeros(tot_pts, 1);
t2_int_vec = zeros(tot_pts, 1);



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
                
                f_local = f_dep_xt(coeff_vec_u,qd_pt_x,qd_pt_t, pt_i, pt_j);
                
                sigma_local = dot(coeff_vec_sigma, bf_eval);
                sigma_x_local = dg_x(coeff_vec_sigma, qd_pts(pt_j));    % to compute sigma_x term
                u_x_local = dg_x(coeff_vec_u, qd_pts(pt_j));            % to compute sigma_x term
                u_t_local = dg_t(coeff_vec_u, qd_pts(pt_i));
                
                % iterate though basis functions for each element

                for bf_j=1:4
                    
                    t2_int_vec(curr_ind(bf_j)) = t2_int_vec(curr_ind(bf_j)) ...
                                            + 0.5*c2*hxht*(diff_const*u_x_local - sigma_local)^2*qd_weights(pt_i)*qd_weights(pt_j);
                                        
                    t1_int_vec(curr_ind(bf_j)) = t1_int_vec(curr_ind(bf_j)) ...
                                            + 0.5*c1*hxht*(u_t_local - sigma_x_local - f_local)^2*qd_weights(pt_i)*qd_weights(pt_j);
               
                end
                   
            end
        end
     end
    
end

% *************** take out boundary points? ******************** %

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

J_eval = 0;

if(strcmp(bdy_cond, 'Dirichlet'))
    J_eval = sum(t1_int_vec(inner_ind_sigma)) + sum(t2_int_vec(inner_ind_u));
end

if(strcmp(bdy_cond, 'Neumann'))
   J_eval = sum(t1_int_vec(inner_ind_sigma)) + sum(t2_int_vec(inner_ind_u));
end




end