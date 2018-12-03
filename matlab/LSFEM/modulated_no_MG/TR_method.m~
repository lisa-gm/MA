% ******************************************************** %
% *************** TRUST REGION METHOD ******************** % 
% ******************************************************** %


function [sigma, ] = TR_method(sigma, u, TR_rad, TR_rad_max, eta, max_iter)


conv_c = max_iter;

eps = 10^(-6);

for k=2:max_iter
    
update_inner = TR_Dogleg(update_inner, grad_J_inner, Hess_J_inner, TR_rad);
update([inner_ind_sigma, tot_pts+inner_ind_u])=update_inner;

update_sigma = update(1:tot_pts);
update_u = update(tot_pts+1:end);

% compute rho
J_eval_old = J_eval_lsfem(hx, Nx_elem, ht, Nt_elem, sigma, u, c1, c2, diff_const);
J_eval_new = J_eval_lsfem(hx, Nx_elem, ht, Nt_elem, sigma+update_sigma, u+update_u, c1, c2, diff_const);

mp = J_eval_old + update_inner'*grad_J_inner + 0.5*update_inner'*Hess_J_inner*update_inner;
rho = (J_eval_old - J_eval_new) / (J_eval_old - mp);
    
if(rho < 0.25)
    TR_rad = 0.25*TR_rad;
elseif(and(rho > 0.75, norm(update_inner) == TR_rad))
    TR_rad = min(2*TR_rad, TR_rad_max);
end
    
if(rho > eta)
    u = u + update_u;
    sigma = sigma + update_sigma;
end


end

% for k=2:max_iter
%     fx = f(x); 
%     
%     g = grad_f(x); 
%     B = hess_f(x); 
%      
%     pk = TR_Dogleg(g, B, TR_rad);
%    
%     % compute rho
%     fp = f(x + pk);
%     mp = fx + pk'*g + 0.5*pk'*B*pk;
%     rho = (fx - fp) / (fx - mp);
%     rho_ratio(k) = rho;
%     
%     if(rho < 0.25)
%         TR_rad = 0.25*TR_rad;
%     elseif(and(rho > 0.75, norm(pk) == TR_rad))
%         TR_rad = min(2*TR_rad, TR_rad_max);
%     end
%     
%     if(rho > eta)
%         x = x + pk;
%     end
%     
%     delta_rad(k) = TR_rad;
%     
%     iter(:, k) = x;
%     
%     if(norm(g) < eps)
%         iter = iter(:, 1:k);
%         conv_c = k;
%         rho_ratio = rho_ratio(1:k);
%         f_coll(k) = f(x);
%         f_coll = f_coll(1:k);
%         delta_rad = delta_rad(1:k);
%        break;
%     end
    
f_coll(k) = f(x);
end

end