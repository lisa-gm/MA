% ******************************************************** %
% ***************** NEWTON STEP ************************** %
% ******************************************************** %


function [sigma, u, norm_grad_J, alpha] = newton_step_lsfem(S, Nx_elem, T, Nt_elem, u0, bdy_cond, bdy_left, bdy_right, sigma, u, J_ss_lin, J_su_lin, J_us_lin, J_uu_lin, c1, c2, diff_const, eps, damp_fact, iter)

%%%%%%%%%%%%%%% SET UP %%%%%%%%%%%%%%
int_space = [0,S];
hx = (int_space(2) - int_space(1))/Nx_elem;
Nx_pts = Nx_elem+1;

int_time = [0,T];
ht = (int_time(2) - int_time(1))/Nt_elem;
Nt_pts = Nt_elem+1;

tot_pts = Nx_pts*Nt_pts;

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

%%%%%%%%%%%%%%%%% computing the gradient %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

grad_J = grad_J_eval(hx, Nx_elem, ht, Nt_elem, bdy_cond, sigma, u, J_ss_lin, J_su_lin, J_us_lin, J_uu_lin, c1);

%%%%%%%%%%%%%%%%%% putting things together: HESSIAN %%%%%%%%%%%%%%%%%%%

Hess_J = Hess_J_eval(hx, Nx_elem, ht, Nt_elem, bdy_cond, sigma, u, J_ss_lin, J_su_lin, J_us_lin, J_uu_lin, c1);


Hess_J_inner = Hess_J([inner_ind_sigma, tot_pts+inner_ind_u], [inner_ind_sigma, tot_pts+inner_ind_u]);
grad_J_inner = grad_J([inner_ind_sigma, tot_pts+inner_ind_u]);

%%%%%%%%%%%%%%%%%%%%%%% SOLVE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
update_inner = - Hess_J_inner \ grad_J_inner;

update = zeros(2*tot_pts, 1);
update([inner_ind_sigma, tot_pts+inner_ind_u])=update_inner;
%%%%%%%%%%%%%% DETERMINE FURTHER UPDATE USING TR %%%%%%%%%%%%%%

max_iter = 10;
eps = 10^(-6);

% choose TR_rad_max
% TR_rad_max = 2;
% TR_rad = 0.75;
% eta = 0.2;
% 
% for k=2:max_iter
%     
% update_inner = TR_Dogleg(update_inner, grad_J_inner, Hess_J_inner, TR_rad);
% update([inner_ind_sigma, tot_pts+inner_ind_u])=update_inner;
% 
% update_sigma = update(1:tot_pts);
% update_u = update(tot_pts+1:end);
% 
% % compute rho
% J_eval_old = J_eval_lsfem(hx, Nx_elem, ht, Nt_elem, sigma, u, c1, c2, diff_const);
% J_eval_new = J_eval_lsfem(hx, Nx_elem, ht, Nt_elem, sigma+update_sigma, u+update_u, c1, c2, diff_const);
% 
% mp = J_eval_old + update_inner'*grad_J_inner + 0.5*update_inner'*Hess_J_inner*update_inner;
% rho = (J_eval_old - J_eval_new) / (J_eval_old - mp);
%     
% if(rho < 0.25)
%     TR_rad = 0.25*TR_rad;
% elseif(and(rho > 0.75, norm(update_inner) == TR_rad))
%     TR_rad = min(2*TR_rad, TR_rad_max);
% end
%     
% if(rho > eta)
%     u = u + update_u;
%     sigma = sigma + update_sigma;
% end
% 
% 
% end

%%%%%%%%%%%%%%%%%%%%%% CHECK UPDATE %%%%%%%%%%%%%%%%%%%%%%%%%%%
temp = dot(update_inner, grad_J_inner);
%fprintf(' sign(dot( update, grad J)) : %d \n', sign(temp));
if(temp > 0)
    % use opposite direction since problem not convex in current direction
    update_inner = - update_inner;
    fprintf('not convex \n');
end

update([inner_ind_sigma, tot_pts+inner_ind_u])=update_inner;

%%%%%%%%%%%%%%%%%%%% LINE SEARCH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha_0 = 1;
c = 0.5;
rho = 0.9;

[alpha, sigma, u, J_eval] = backtracking(alpha_0, c, rho, hx, Nx_elem, ht, Nt_elem, bdy_cond, sigma, u, update, J_ss_lin, J_su_lin, J_us_lin, J_uu_lin, c1, c2, diff_const);

if(mod(iter, 10) == 1)
    J_eval_vec = plot_J(hx, Nx_elem, ht, Nt_elem, sigma, u, c1, c2, diff_const, iter);
end

% compute J of current solution to compare with potential updates
% J_eval_old = J_eval_lsfem(S, Nx_elem, T, Nt_elem, sigma, u, c1, c2, diff_const);
% 
% inner_it = 0;

% u_temp = u;
% sigma_temp = sigma;

% 
% while(inner_it < 60)    
%     
%     sigma_temp(inner_ind_sigma) = sigma(inner_ind_sigma)...
%             + damp_fact*update_inner(1:length(inner_ind_sigma));
%     u_temp(inner_ind_u) = u(inner_ind_u) ...
%             + damp_fact*update_inner(length(inner_ind_sigma)+1:end);
%     
%     J_eval_temp = J_eval_lsfem(S, Nx_elem, T, Nt_elem, sigma_temp, u_temp, c1, c2, diff_const);
% 
%     fprintf('damping factor = %g, new J = %g, old J = %g\n', damp_fact,J_eval_temp,J_eval_old);
%     
%     if J_eval_temp < J_eval_old
%       break;
%     end
%     
%     damp_fact = damp_fact/2;
%     inner_it = inner_it + 1;
%     if(inner_it == 50)
%        fprintf('max iter reached!');
%     end
% end

norm_grad_J = norm(grad_J_eval(hx, Nx_elem, ht, Nt_elem, bdy_cond, sigma, u, J_ss_lin, J_su_lin, J_us_lin, J_uu_lin, c1));

fprintf('iter: %g, norm J = %1.15d, norm grad J = %1.15d, norm update = %g, alpha = %g\n', iter, J_eval, norm_grad_J, norm(alpha*update), alpha);

% enforcing bdy conditions on sigma & u

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
