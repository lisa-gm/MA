%%%%%% V cycle %%%%%%

% restrict ourselves to grid sizes of type 4,8, 16, 32 , ... for now

% J gives the number of coarsening levels

%%%%%% V cycle %%%%%%

% restrict ourselves to grid sizes of type 4,8, 16, 32 , ... for now

% J gives the number of coarsening levels
function u = V_cycle(A, f, u0, J, max_iter, smoother)

N = length(f);
eps = 10^(-12);

max_iter_sm = 3;

% set up matrices for MG
[L_h, R, I, u0, rhs] = set_up_mg(A, f, J);

    % check  if N/(2^J) < 1 then error
if (N/(2^J) < 1)
    fprintf('too many levels');
    return;
end

% write function for smoothening loop
for k=1:max_iter
    
    % overall loop, going up coarsening steps
for j=1:J 
  [rhs{j+1}, u0{j}] = step_fine_to_coarse(L_h{j}, R{j}, rhs{j}, u0{j}, max_iter_sm, smoother);
end

    % use direct solver
    u0{J+1} = L_h{J+1}\(rhs{J+1});
    
    % walk back from coarse to fine adding in corrections
    % for J u0{J+1} = e{J+1}
    % but then 
    % u0{j} = u0{j} + e{j}
    % where e{j} is interpolated correction from previous step
    
    for j=J:-1:1
        u0{j} = u0{j} + I{j}*u0{j+1}; 
    end
    
    
% do some steps of post smoothening
    if(strcmp(smoother, 'Jacobi'))
        [u0{1}, ~] = JacobiSolve(A, f, u0{1}, max_iter_sm);
    end
    
    if(strcmp(smoother, 'GaussSeidel'))
        [u0{1}, ~] = GaussSeidelSolve(A, f, u0{1}, max_iter_sm);
    end
    
    if(strcmp(smoother, 'ConjugateGradient'))
        [u0{1}, ~] = CGSolve(A, f, u0{1}, max_iter_sm);
        
    end
    
    
% u0{1} should give final solution of V cycle
err = norm(f-A*u0{1});

if(err < eps)
    %fprintf(' multigrid converged after %d V-cycles\n', k);
    break;
end

end
u = u0{1};
end