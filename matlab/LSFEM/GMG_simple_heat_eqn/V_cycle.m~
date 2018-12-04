% ******************************************************************* %
% ******************** V - CYCLE GEOMEMTRIC MG ********************** %
% ******************************************************************* %


% J gives the number of coarsening levels
function u = V_cycle(A, f, bdy_ind, J, max_iter, smoother, eps)

N = length(f);
%eps = 10^(-9);

max_iter_sm = 2;

% set up matrices for MG
[L_h, R, I, u0, rhs] = set_up_mg(A, f, J, bdy_ind);
res_after_pre_sm = {};

    % check  if N/(2^J) < 1 then error
if (N/(2^J) < 1)
    fprintf('too many levels');
    return;
end


% write function for smoothening loop
for k=1:max_iter
    
    % overall loop, going up coarsening steps
for j=1:J 
  %fprintf('on level: %d ', j);
  [rhs{j+1}, u0{j}] = step_fine_to_coarse(L_h{j}, R{j}, rhs{j}, u0{j}, max_iter_sm, smoother, j);
  res_after_pre_sm{j} = rhs{j} - L_h{j}*u0{j};
end

    % use direct solver
    u0{J+1} = L_h{J+1}\(rhs{J+1});
    
    % walk back from coarse to fine adding in corrections
    % for J u0{J+1} = e{J+1}
    % but then 
    % u0{j} = u0{j} + e{j}
    % where e{j} is interpolated correction from previous step
    
    for j=J:-1:1
        % alpha scaling parameter 
        c = I{j}*u0{j+1};
        alpha = 1;
        %alpha = dot(c, res_after_pre_sm{j})/dot(L_h{j}*c, c);
        u0{j} = u0{j} + alpha*c;
        
        if(j == 1)
        fprintf('energy before post-smoothing : %d\n', abs(u0{j}'*L_h{j}*u0{j}-rhs{j}'*u0{j}));
        end
        
        u0{j} = smoothing(L_h{j}, rhs{j}, u0{j}, max_iter_sm, smoother);
        
        if(j == 1)
        fprintf('energy after post-smoothing : %d\n', abs(u0{j}'*L_h{j}*u0{j}-rhs{j}'*u0{j}));
        end
        
    end    
    
% u0{1} should give final solution of V cycle
err = norm(f-A*u0{1});
%fprintf('iter : %d, norm residual : %d\n', k, err);

if(err < eps)
    fprintf(' multigrid converged after %d V-cycles\n', k);
    break;
end

end
u = u0{1};
end