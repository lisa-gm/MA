%%%%% helper function to go from fine to coarse %%%%

function [rhs, u0] = step_fine_to_coarse(L_h, R, rhs, u0, max_iter_sm, smoother, level)
%%%%%%%% loop to do pre smoothening


    %%% STEP 1: smooth L_h * u_approx = f by applying 
    % iterative scheme on fine grid
    if( level == 1)
     fprintf('norm res before pre-smoothing : %d\n', norm(L_h*u0 - rhs));
    end
    
    u0 = smoothing(L_h, rhs, u0, max_iter_sm, smoother);

    if( level == 1)
     fprintf('norm res after pre-smoothing : %d\n', norm(L_h*u0 - rhs));
    end
    
    %%%%% STEP 2: compute residual, then restrict residual
    % to coarse grid
    r = rhs - L_h*u0;
    
    % now restrict r{j}, and hence define new RHS
    rhs = R*r;    
end