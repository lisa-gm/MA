% ************************************************************ %
% *************** BACKTRACKING LINE SEARCH ******************* %
% ************************************************************ %

function [alpha, sigma, u, J_temp] = backtracking(alpha_0, c, rho, hx, Nx_elem, ht, Nt_elem, bdy_cond, sigma, u, update, J_ss_lin, J_su_lin, J_us_lin, J_uu_lin, c1, c2, diff_const)
    alpha = alpha_0;
    
    update_sigma = update(1:length(update)/2);
    update_u = update(length(update)/2+1:end);
    
    [J_old, ~] = J_eval_lsfem(hx, Nx_elem, ht, Nt_elem, sigma, u, c1, c2, diff_const);
    grad_J_old = grad_J_eval(hx, Nx_elem, ht, Nt_elem, bdy_cond, sigma, u, J_ss_lin, J_su_lin, J_us_lin, J_uu_lin, c1);

    [J_temp, ~] = J_eval_lsfem(hx, Nx_elem, ht, Nt_elem, sigma+alpha*update_sigma, u+alpha*update_u, c1, c2, diff_const);


    while( J_temp > (J_old + c*alpha*(grad_J_old)'*update))
        alpha = alpha * rho;
        [J_temp, ~] = J_eval_lsfem(hx, Nx_elem, ht, Nt_elem, sigma+alpha*update_sigma, u+alpha*update_u, c1, c2, diff_const);
    end
    
    norm(alpha*[update_sigma; update_u]);
    
    sigma = sigma + alpha*update_sigma;
    u = u + alpha*update_u;
end