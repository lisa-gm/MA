%%%%% helper function to go from fine to coarse %%%%

function [rhs, u0] = step_fine_to_coarse(L_h, R, rhs, u0, max_iter_sm, smoother)
%%%%%%%% loop to do pre smoothening


    %%% STEP 1: smooth L_h * u_approx = f by applying 
    % iterative scheme on fine grid
    
    u0 = smoothing(L_h, rhs, u0, max_iter_sm, smoother);

    % smoother_list = ['Jacobi', 'GaussSeidel', 'ConjugateGradient', 'GMRES'];

%     if(ismember(smoother, smoother_list) == 0)
%         disp(' Invalid Smoother! ')
%         return;
%     end
%     
%     if(strcmp(smoother, 'Jacobi'))
%         [u0, ~] = JacobiSolve(L_h, rhs, u0, max_iter_sm);
%     end
%     
%     if(strcmp(smoother, 'GaussSeidel'))
%         u0 = GaussSeidelSolve(L_h, rhs, u0, max_iter_sm);
%     end
%     
%     if(strcmp(smoother, 'ConjugateGradient'))
%         [u0, ~] = CGSolve(L_h, rhs, u0, max_iter_sm);
%     end
%     
%     if(strcmp(smoother, 'GMRES'))
%         
%     end

    %%%%% STEP 2: compute residual, then restrict residual
    % to coarse grid
    r = rhs - L_h*u0;
    
    % now restrict r{j}, and hence define new RHS
    rhs = R*r;
end