function [L_h, R, I, u0, rhs] = set_up_mg(Nx_elem_list, hx_list, Nt_elem_list, ht_list, c1, c2, diff_const, bdy_cond, bdy_left, bdy_right, levels);
S = hx_list(1)*Nx_elem_list(1);
T = ht_list(1)*Nt_elem_list(1);

% interpolation, restriction, operators, initial guesses
I = {};
R = {};

L_h = {};
u0 = {};
rhs = {};

% set up matrices for grid transfer 
% only choose every other entry
    for j=1:J        
    sub_size = 2*(Nx_elem_list(j)+1)*(Nt_elem_list(j)+1);
    
    % setting interpolation & restriction operator
    temp = ones(sub_size,1);
    mat_temp = spdiags([temp, 2*temp, temp], -1:1, sub_size, sub_size);
    if(j == 1)
        I{j} = (1/2*mat_temp(2:2:end, :))';
        I{j}(bdy_ind, :) = 0;
        %R{j} = 1/4*mat_temp(2:2:end, :);
        R{j} = I{j}';
    else
        I{j} = (1/2*mat_temp(2:2:end, :))';
        %R{j} = 1/4*mat_temp(2:2:end, :);
        R{j} = I{j}';
   end
    
    % to store intermediate results
    u0{j} = zeros(sub_size,1);  
            
    % reconstruct original operator just coarser grid    
     
    [J_ss_lin, J_su_lin, J_us_lin, J_uu_lin] = lin_part_mat_no_bdy_lsfem (S, Nx_elem, T, Nt_elem, c1 , c2, diff_const, 'none');
    L_h{j} = [J_ss_lin, J_su_lin; J_us_lin, J_uu_lin];
    
    rhs{j} = zeros(sub_size, 1);
    end
    
    end

    end