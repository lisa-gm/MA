% ******************************************************************* %
% ******************** V - CYCLE GEOMEMTRIC MG ********************** %
% ******************************************************************* %


% J gives the number of coarsening levels
function [sol, it_sol_mg, it_res_mg] = V_cycle_galerkin_assembly_P(Nx_elem_list, hx_list, Nt_elem_list, ht_list, Hess_J, f, sol_init, bdy_cond, u0, bdy_left, bdy_right, levels, max_iter, smoother, eps)
max_iter_sm = 3;

% up until now only works for 2 level system, therefore 
if(levels > 2)
    fprintf('ATM ONLY 2 LEVELS POSSIBLE ... SOZ');
    pause;
end


% get interpolation matrices
[I, R] = set_up_interpolation_op_SP_TIME(Nx_elem_list, Nt_elem_list, bdy_cond, levels);

rhs = {};
L_h = {};
u = {};
    
bdy_ind_list = {};

% overall loop, going up coarsening steps
for l=1:levels   
  % get Hessian on each level
  if(l == 1)
    L_h{l} = Hess_J;
    rhs{l} = f;
    u{l} = sol_init;
    
    bdy_ind_list{l} = get_bdy_ind(Nx_elem_list(l)+1, Nt_elem_list(l)+1, bdy_cond);
   [rhs{l}, L_h{l}] = apply_bdy_cond(Nx_elem_list(l)+1, Nt_elem_list(l)+1, L_h{l}, rhs{l}, bdy_cond, u0, bdy_left, bdy_right, l);
  
  else
      L_h{l} = R{l-1}*L_h{l-1}*I{l-1};
      rhs{l} = zeros(size(L_h{l},1),1);

       
   [rhs{l}, L_h{l}] = apply_bdy_cond(Nx_elem_list(l)+1, Nt_elem_list(l)+1, L_h{l}, rhs{l}, bdy_cond, u0, bdy_left, bdy_right, l);
    bdy_ind_list{l} = get_bdy_ind(Nx_elem_list(l)+1, Nt_elem_list(l)+1, bdy_cond);
  end
  
end

% construct smoother once here, then never again. pass through.
block_size_s = 2;
block_size_t = 2;
P = construct_block_Jac_mat(Nx_elem_list(l)+1, Nt_elem_list(l)+1, L_h{1}, block_size_s, block_size_t);

it_sol_mg = zeros(length(sol_init), max_iter+1);
it_sol_mg(:,1) = sol_init;

it_res_mg = zeros(1, max_iter+1);
it_res_mg(1,1) = norm(L_h{1}*sol_init - f);

% write function for smoothening loop
for k=1:max_iter
   
  for l=1:levels-1 

      [rhs{l+1}, u{l}] = step_fine_to_coarse(Nx_elem_list(l)+1, Nt_elem_list(l)+1, L_h{l}, R{l}, rhs{l}, u{l}, max_iter_sm, smoother, l); 

      % remove bdy conditions from residual
      rhs{l+1}(bdy_ind_list{l+1}) = 0;
  end
  
  % use direct solver
  u{levels} = L_h{levels}\(rhs{levels});
    
    for l=levels-1:-1:1
        % alpha scaling parameter 
        c = I{l}*u{l+1};
        
        % set c to zero on bdy points
        c(bdy_ind_list{l}) = 0;

        u{l} = u{l} + c;
        % fprintf('norm res before post-smoothing : %d\n', norm(L_h{l}*u{l} - rhs{l}));

        u{l} = smoothing(Nx_elem_list(l)+1, Nt_elem_list(l)+1 ,L_h{l}, rhs{l}, u{l}, max_iter_sm, smoother, P);
        
        if(l == 1)
        fprintf('norm res after post-smoothing : %d\n', norm(L_h{l}*u{l} - rhs{l}));
        end
        
    end    
    
it_sol_mg(:,k+1) = u{1};

err = norm(L_h{1}*u{1} - rhs{1});
it_res_mg(1,k+1) = err;


sol = u{1};


    if(err < eps)
        fprintf('\n');
        fprintf('norm residual : %d\n', err);
        fprintf(' multigrid P assembly converged after %d V-cycles\n', k);
        it_sol_mg = it_sol_mg(:,1:k+1);
        it_res_mg = it_res_mg(1,1:k+1);
        break;
    end
end  % end max iter

end