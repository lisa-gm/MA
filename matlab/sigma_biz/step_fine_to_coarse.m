%%%%% helper function to go from fine to coarse %%%%

function [rhs, u0] = step_fine_to_coarse(L_h, R, rhs, u0, max_iter_sm, smoother)
%%%%%%%% loop to do pre smoothening
size(L_h);
size(R);
size(rhs);
size(u0);

    %%% STEP 1: smooth L_h * u_approx = f by applying 
    % iterative scheme on fine grid
    smoother_list = ['Jacobi', 'GaussSeidel', 'ConjugateGradient'];

    if(ismember(smoother, smoother_list) == 0)
        disp(' Invalid Smoother! ')
        return;
    end
    
    if(strcmp(smoother, 'Jacobi'))
        [u0, ~] = JacobiSolve(L_h, rhs, u0, max_iter_sm);
    end
    
    if(strcmp(smoother, 'GaussSeidel'))
        [u0, ~] = GaussSeidelSolve(L_h, rhs, u0, max_iter_sm);
    end
    
    if(strcmp(smoother, 'ConjugateGradient'))
        [u0, ~] = CGSolve_2(L_h, rhs, u0, max_iter_sm);
    end

    %%%%% STEP 2: compute residual, then restrict residual
    % to coarse grid
    r = rhs - L_h*u0;
    
    % now restrict r{j}, and hence define new RHS
    rhs = R*r;
end