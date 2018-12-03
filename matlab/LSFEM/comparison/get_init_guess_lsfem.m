% ******************************************************** %
% *************** GET BETTER INITIAL GUESS *************** %
% ******************************************************** %

function [sigma_init, u_init] = get_init_guess_lsfem(S, Nx_elem, T, Nt_elem, u0, bdy_cond, bdy_left, bdy_right, a0, a1, J_ss_lin, J_su_lin, J_us_lin, J_uu_lin, c1, c2, diff_const)

int_space = [0,S];
hx = (int_space(2) - int_space(1))/Nx_elem;
Nx_pts = Nx_elem+1;
x_vec = linspace(int_space(1), int_space(2), Nx_pts);

int_time = [0,T];
ht = (int_time(2) - int_time(1))/Nt_elem;
Nt_pts = Nt_elem+1;
t_vec = linspace(int_time(1), int_time(2), Nt_pts);

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% CONSRUCTING THE OPERATORS %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
                    J_s_a0(curr_ind(bf_j)) = J_s_a0(curr_ind(bf_j)) - c1*hxht*grad_tau_x_eval(bf_j)*a0(qd_pt_x, qd_pt_t)*qd_weights(pt_i)*qd_weights(pt_j);            %assemble rhs                                                  %rhs
                    J_u_a0(curr_ind(bf_j)) = J_u_a0(curr_ind(bf_j)) + c1*hxht*grad_v_t_eval(bf_j)*a0(qd_pt_x, qd_pt_t)*qd_weights(pt_i)*qd_weights(pt_j);
                    J_u_a0(curr_ind(bf_j)) = J_u_a0(curr_ind(bf_j)) - c1*hxht*bf_eval(bf_j)*a1(qd_pt_x, qd_pt_t)*a0(qd_pt_x, qd_pt_t)*qd_weights(pt_i)*qd_weights(pt_j);
                 
                    for bf_i=1:4
                      
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
J_us_init(1:Nx_pts,:) = 0;

J_uu_init(1:Nx_pts, 1:Nx_pts) = eye(Nx_pts, Nx_pts);

if(strcmp(bdy_cond, 'Dirichlet'))
    J_uu_init(1:Nx_pts:end) = 0;
    J_uu_init(Nx_pts:Nx_pts:end) = 0;

    J_us_init(1:Nx_pts:end) = 0;
    J_us_init(Nx_pts:Nx_pts:end) = 0;

    J_uu_init(1:Nx_pts:end, 1:Nx_pts:end) = eye(Nt_pts, Nt_pts);
    J_uu_init(Nx_pts:Nx_pts:end, Nx_pts:Nx_pts:end) = eye(Nt_pts, Nt_pts);
end

if(strcmp(bdy_cond, 'Neumann'))
    J_ss_init(1:Nx_pts:end) = 0;
    J_ss_init(Nx_pts:Nx_pts:end) = 0;

    J_su_init(1:Nx_pts:end) = 0;
    J_su_init(Nx_pts:Nx_pts:end) = 0;

    J_ss_init(1:Nx_pts:end, 1:Nx_pts:end) = eye(Nt_pts, Nt_pts);
    J_ss_init(Nx_pts:Nx_pts:end, Nx_pts:Nx_pts:end) = eye(Nt_pts, Nt_pts);
end

% ASSEMBLE MATRIX 
J_init_guess = [J_ss_init, J_su_init; J_us_init, J_uu_init];
block_size = size(J_ss_init,1);

% RHS 
F_init = [J_s_a0; J_u_a0];
F_init(block_size+1:block_size+Nx_pts) = u0;

if(strcmp(bdy_cond, 'Dirichlet'))
    F_init(block_size+1:Nx_pts:end) = bdy_left;
    F_init(block_size+Nx_pts:Nx_pts:end) = bdy_right;
end

if(strcmp(bdy_cond, 'Neumann'))
    F_init(1:Nx_pts:block_size) = -bdy_left;
    F_init(Nx_pts:Nx_pts:block_size) = bdy_right;
end

%%%%%%%%%%%%%%%%%%%%% SOLVE  %%%%%%%%%%%%%%%%%%%%%%
sol = J_init_guess \ F_init;

% INITIAL GUESS
sigma = sol(1:block_size);
u = sol(block_size+1:end);

% reshape
u_mat_init = zeros(Nx_pts, Nt_pts);
sigma_mat_init = zeros(Nx_pts, Nt_pts);


for ts=1:Nt_pts
    u_mat_init(:, ts) = u((ts-1)*Nx_pts+1:ts*Nx_pts);
    sigma_mat_init(:, ts) = sigma((ts-1)*Nx_pts+1:ts*Nx_pts);
    
end

% figure
% mesh(x_vec, t_vec(1:end), transpose(u_mat_init(:,1:end)));
% xlabel('space');
% ylabel('time');
% zlabel('u');
% %zlim([min(u), max(u)]);
% title('heat equation LSFEM w/ lin part of f');
% 
% figure
% mesh(x_vec, t_vec(1:end), transpose(sigma_mat_init(:,1:end)));
% xlabel('space');
% ylabel('time');
% zlabel('u');
% %zlim([min(u), max(u)]);
% title('SIGMA heat equation LSFEM w/ lin part of f');

%sigma = reshape(sigma_exact, [tot_pts, 1]);
%u = reshape(u_exact, [tot_pts, 1]);

if(strcmp(bdy_cond, 'Neumann'))
    sigma(1:Nx_pts:end) = -bdy_left;
    sigma(Nx_pts:Nx_pts:end) = bdy_right;
end

sigma_init = sigma;
u_init = u;

end