% ******************************************************************* %
% ******************** V - CYCLE GEOMEMTRIC MG ********************** %
% ******************************************************************* %


% J gives the number of coarsening levels
function [sol, it_sol_mg, it_res_mg] = V_cycle(Nx_elem_list, hx_list, Nt_elem_list, ht_list, f, sol_init, c1, c2, diff_const, bdy_cond, u0, bdy_left, bdy_right, levels, max_iter, smoother, eps)
S = Nx_elem_list(1)*hx_list(1);
T =  Nt_elem_list(1)*ht_list(1);

max_iter_sm = 2;

% get interpolation matrices
%[I, R, Nx_elem_list, Nt_elem_list] = set_up_interpolation_op_SP(Nx_elem_list, Nt_elem_list, bdy_cond, levels);
[I, R] = set_up_interpolation_op_SP_TIME(Nx_elem_list, Nt_elem_list, bdy_cond, levels);


Nx_pts_list = Nx_elem_list +1;
Nt_pts_list = Nt_elem_list +1;

rhs = {};
L_h = {};
u = {};
    
bdy_ind_list = {};

% overall loop, going up coarsening steps
for l=1:levels   
  % get Hessian on each level
  if(l == 1)
    [J_ss_lin, J_su_lin, J_us_lin, J_uu_lin] = lin_part_mat_no_bdy_lsfem (S, Nx_elem_list(l), T, Nt_elem_list(l), c1 , c2, diff_const);
    L_h{l} = [J_ss_lin, J_su_lin; J_us_lin, J_uu_lin];
    rhs{l} = f;
    u{l} = sol_init;
    
    bdy_ind_list{l} = get_bdy_ind(Nx_elem_list(l)+1, Nt_elem_list(l)+1, bdy_cond);
   [rhs{l}, L_h{l}] = apply_bdy_cond(Nx_elem_list(l)+1, Nt_elem_list(l)+1, L_h{l}, rhs{l}, bdy_cond, u0, bdy_left, bdy_right, l);
  
  else
%     L_h{l} = R{l-1}*L_h{l-1}*I{l-1}; 
%     rhs{l} = zeros(size(L_h{l},1),1);
%     u{l} = zeros(size(L_h{l},1),1);
%     [~, L_h{l}] = apply_bdy_cond(Nx_elem_list(l)+1, Nt_elem_list(l)+1, L_h{l}, rhs{l}, bdy_cond, u0, bdy_left, bdy_right, l);  

    [J_ss_lin, J_su_lin, J_us_lin, J_uu_lin] = lin_part_mat_no_bdy_lsfem (S, Nx_elem_list(l), T, Nt_elem_list(l), c1 , c2, diff_const);
    L_h{l} = [J_ss_lin, J_su_lin; J_us_lin, J_uu_lin];
    rhs{l} = zeros(size(L_h{l},1),1);
    u{l} = zeros(size(L_h{l},1),1);
    
    
   [rhs{l}, L_h{l}] = apply_bdy_cond(Nx_elem_list(l)+1, Nt_elem_list(l)+1, L_h{l}, rhs{l}, bdy_cond, u0, bdy_left, bdy_right, l);
    bdy_ind_list{l} = get_bdy_ind(Nx_elem_list(l)+1, Nt_elem_list(l)+1, bdy_cond);
  end
  
end

it_sol_mg = zeros(length(sol_init), max_iter+1);
it_sol_mg(:,1) = sol_init;

it_res_mg = zeros(1, max_iter+1);
it_res_mg(1,1) = norm(L_h{1}*sol_init - f);

% write function for smoothening loop
for k=1:max_iter
    
    % reshape
%     u_temp = u{1}(Nx_pts_list(1)*Nt_pts_list(1)+1:end);
%     u_mat_temp = zeros(Nx_pts_list(1), Nt_pts_list(1));
% 
%     for ts=1:Nt_pts_list(1)
%         u_mat_temp(:, ts) = u_temp((ts-1)*Nx_pts_list(1)+1:ts*Nx_pts_list(1));
%     end
% 
%     figure(1);
%     mesh(0:Nx_elem_list(1), 0:Nt_elem_list(1), transpose(u_mat_temp));
%     xlabel('space');
%     ylabel('time');
%     zlabel('u temp');
%     %zlim([min(u), max(u)]);
%     title('current sol');
   
  for l=1:levels-1 
  %fprintf('on level: %d ', j);
  [rhs{l+1}, u{l}] = step_fine_to_coarse(L_h{l}, R{l}, rhs{l}, u{l}, max_iter_sm, smoother, l); 
  % now remove boundary conditions
      % reshape
      
%     if(l == 1)
%         u_temp = u{1}(Nx_pts_list(1)*Nt_pts_list(1)+1:end);
%         u_mat_temp = zeros(Nx_pts_list(1), Nt_pts_list(1));
% 
%         for ts=1:Nt_pts_list(1)
%             u_mat_temp(:, ts) = u_temp((ts-1)*Nx_pts_list(1)+1:ts*Nx_pts_list(1));
%         end
% 
%         figure(2);
%         mesh(0:Nx_elem_list(1), 0:Nt_elem_list(1), transpose(u_mat_temp));
%         xlabel('space');
%         ylabel('time');
%         zlabel('u temp');
%         %zlim([min(u), max(u)]);
%         title('after pre smoothing');
%     end
    
        if(l == 1)
        sigma_temp = u{1}(1:Nx_pts_list(1)*Nt_pts_list(1));
        sigma_mat_temp = zeros(Nx_pts_list(1), Nt_pts_list(1));

        for ts=1:Nt_pts_list(1)
            sigma_mat_temp(:, ts) = sigma_temp((ts-1)*Nx_pts_list(1)+1:ts*Nx_pts_list(1));
        end

        figure(10);
        mesh(0:Nx_elem_list(1), 0:Nt_elem_list(1), transpose(sigma_mat_temp));
        xlabel('space');
        ylabel('time');
        zlabel('sigma temp');
        %zlim([min(u), max(u)]);
        title('after pre smoothing');
        end
 
  
  rhs{l+1}(bdy_ind_list{l+1}) = 0;
  end

  % use direct solver
  u{levels} = L_h{levels}\(rhs{levels});
  
  u_vec_c = u{levels}(Nx_pts_list(end)*Nt_pts_list(end)+1:end);
  u_mat_c = zeros(Nx_pts_list(end), Nt_pts_list(end));

    for ts=1:Nt_pts_list(end)
        u_mat_c(:, ts) = u_vec_c((ts-1)*Nx_pts_list(end)+1:ts*Nx_pts_list(end));
    end

    figure(3);
    mesh(1:Nx_pts_list(end), 1:Nt_pts_list(end), transpose(u_mat_c));
    xlabel('space');
    ylabel('time');
    zlabel('u');
    %zlim([min(u), max(u)]);
    title('coarse grid correction');
%     
%     sigma_vec_c = u{levels}(1:Nx_pts_list(end)*Nt_pts_list(end));
%     sigma_mat_c = zeros(Nx_pts_list(end), Nt_pts_list(end));
% 
%     for ts=1:Nt_pts_list(end)
%         sigma_mat_c(:, ts) = sigma_vec_c((ts-1)*Nx_pts_list(end)+1:ts*Nx_pts_list(end));
%     end
% 
%     figure(13);
%     mesh(1:Nx_pts_list(end), 1:Nt_pts_list(end), transpose(sigma_mat_c));
%     xlabel('space');
%     ylabel('time');
%     zlabel('u');
%     %zlim([min(u), max(u)]);
%     title('coarse grid correction');
    

  
  
  % walk back from coarse to fine adding in corrections
  % for J u0{J+1} = e{J+1}
  % but then 
  % u0{j} = u0{j} + e{j}
  % where e{j} is interpolated correction from previous step
    
    for l=levels-1:-1:1
        % alpha scaling parameter 
        c = I{l}*u{l+1};
        
        % set c to zero on bdy points
        c(bdy_ind_list{l}) = 0;
        
        if(l==1)
                % reshape
                
            sigma_temp = u{1}(1:Nx_pts_list(1)*Nt_pts_list(1));
            sigma_mat_temp = zeros(Nx_pts_list(1), Nt_pts_list(1));
                       

            c_sigma_temp = c(1:Nx_pts_list(1)*Nt_pts_list(1));
            c_sigma_mat_temp = zeros(Nx_pts_list(1), Nt_pts_list(1));

            for ts=1:Nt_pts_list(1)
                sigma_mat_temp(:, ts) = sigma_temp((ts-1)*Nx_pts_list(1)+1:ts*Nx_pts_list(1));
                c_sigma_mat_temp(:, ts) = c_sigma_temp((ts-1)*Nx_pts_list(1)+1:ts*Nx_pts_list(1));
            end


            figure(4);
            mesh(0:Nx_elem_list(1), 0:Nt_elem_list(1), transpose(c_sigma_mat_temp));
            xlabel('space');
            ylabel('time');
            zlabel('sigma temp');
            %zlim([min(u), max(u)]);
            title('fine grid correction');
            
            figure(5);
            mesh(0:Nx_elem_list(1), 0:Nt_elem_list(1), transpose(sigma_mat_temp));
            xlabel('space');
            ylabel('time');
            zlabel('u temp');
            %zlim([min(u), max(u)]);
            title('sigma before update: sigma mat temp');
            
            figure(6);
            mesh(0:Nx_elem_list(1), 0:Nt_elem_list(1), transpose(sigma_mat_temp+c_sigma_mat_temp));
            xlabel('space');
            ylabel('time');
            zlabel('u temp');
            %zlim([min(u), max(u)]);
            title('current sol: sigma mat temp + correction');
            
        
        end
        
        alpha = 1;
        %alpha = dot(c, res_after_pre_sm{l})/dot(L_h{l}*c, c);
        u{l} = u{l} + alpha*c;
        
        if(l == 1)
        res = L_h{l}*u{l} - rhs{l};
        tot_pts = length(res)/2;
        fprintf('norm res before post-smoothing : %d\n', norm(res));
        %fprintf('norm res sigma before post-smoothing : %d\n', norm(res(1:tot_pts)));
        %fprintf('norm res u before post-smoothing : %d\n', norm(res(tot_pts+1:end)));

        end
        
        u{l} = smoothing(L_h{l}, rhs{l}, u{l}, max_iter_sm, smoother);
        
        if(l == 1)
        fprintf('norm res after post-smoothing : %d\n', norm(L_h{l}*u{l} - rhs{l}));
        end
        
    end    
    
it_sol_mg(:,k+1) = u{1};

err = norm(L_h{1}*u{1} - rhs{1});
it_res_mg(1,k+1) = err;


if(err < eps)
    fprintf('\n');
    fprintf('norm residual : %d\n', err);
    fprintf(' multigrid converged after %d V-cycles\n', k);
    it_sol_mg = it_sol_mg(:,1:k+1);
    it_res_mg = it_res_mg(1,1:k+1);
    break;
end

end

sol = u{1};
end