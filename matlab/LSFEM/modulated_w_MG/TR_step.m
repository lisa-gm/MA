% ******************************************************** %
% *************** TRUST REGION METHOD ******************** % 
% ******************************************************** %


function [sigma, u, norm_grad_J, rho, TR_rad] = TR_step(TR_rad, TR_rad_max, eta, max_iter_TR, S, Nx_elem, T, Nt_elem, u0, bdy_cond, bdy_left, bdy_right, sigma, u, update, grad_J, Hess_J, J_ss_lin, J_su_lin, J_us_lin, J_uu_lin, c1, c2, diff_const, iter)
eps = 10^(-9);
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

% grad_J = grad_J_eval(hx, Nx_elem, ht, Nt_elem, bdy_cond, sigma, u, J_ss_lin, J_su_lin, J_us_lin, J_uu_lin, c1);
% 
% %%%%%%%%%%%%%%%%%% putting things together: HESSIAN %%%%%%%%%%%%%%%%%%%
% 
% Hess_J = Hess_J_eval(hx, Nx_elem, ht, Nt_elem, bdy_cond, sigma, u, J_ss_lin, J_su_lin, J_us_lin, J_uu_lin, c1);
% 
% 
Hess_J_inner = Hess_J([inner_ind_sigma, tot_pts+inner_ind_u], [inner_ind_sigma, tot_pts+inner_ind_u]);
grad_J_inner = grad_J([inner_ind_sigma, tot_pts+inner_ind_u]);
update_inner = update([inner_ind_sigma, tot_pts+inner_ind_u]);

%%%%%%%%%%%%%%%%%%%%%%% SOLVE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
counter = 1;
while true
    
    [J_eval_old, J_eval_old_vec] = J_eval_lsfem(hx, Nx_elem, ht, Nt_elem, sigma, u, c1, c2, diff_const);

    % find best dogleg update without solving LSE again
    update_inner = TR_Dogleg(grad_J_inner, Hess_J_inner, update_inner, TR_rad);
    update([inner_ind_sigma, tot_pts+inner_ind_u])=update_inner;

    update_sigma = update(1:tot_pts,1);
    update_u = update(tot_pts+1:end,1);

    [J_eval_new, J_eval_new_vec] = J_eval_lsfem(hx, Nx_elem, ht, Nt_elem, sigma+update_sigma, u+update_u, c1, c2, diff_const);

    mp = J_eval_old + update_inner'*grad_J_inner + 0.5*update_inner'*Hess_J_inner*update_inner;
    rho = (J_eval_old - J_eval_new) / (J_eval_old - mp);

    if(rho < 0.25)
        TR_rad = 0.25*TR_rad;
    elseif(rho > 0.75)
        TR_rad = min(2*TR_rad, 9e9);
    end

    if(rho > eta)
        u = u + update_u;
        sigma = sigma + update_sigma;
        
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
        
        Hess_J = Hess_J_eval(hx, Nx_elem, ht, Nt_elem, bdy_cond, sigma, u, J_ss_lin, J_su_lin, J_us_lin, J_uu_lin, c1);
        grad_J = grad_J_eval(hx, Nx_elem, ht, Nt_elem, bdy_cond, sigma, u, J_ss_lin, J_su_lin, J_us_lin, J_uu_lin, c1);
        
        Hess_J_inner = Hess_J([inner_ind_sigma, tot_pts+inner_ind_u], [inner_ind_sigma, tot_pts+inner_ind_u]);
        grad_J_inner = grad_J([inner_ind_sigma, tot_pts+inner_ind_u]);
        
        %%%%%%%%%%%%% CALL MULTIGRID FROM HERE %%%%%%%%%%%%%%%%%%%%%%%%%%
        update_inner = - Hess_J_inner \ grad_J_inner;
        

        norm_grad_J = norm(grad_J);
        J_eval = J_eval_lsfem(hx, Nx_elem, ht, Nt_elem, sigma, u, c1, c2, diff_const);
   
        fprintf('iter: %g, norm J = %1.15d, norm grad J = %1.15d, norm update = %g, rho = %g, TR_rad =%g\n', counter, J_eval, norm_grad_J, norm([update_sigma; update_u]), rho, TR_rad);
        counter = counter + 1;
        
        % break out of TR_step to recompute a new Hessian, gradient, update
        %return;
        
        %if(norm_grad_J < 10^(-9))
            %fprintf('convergence criterion reached after %d iterations\n', iter);
            %break;
        %end
        
%         if(counter == 5)
%             break;
%         end
        %break;
    end
        
end



end
