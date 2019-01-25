%%%%% helper function to go from fine to coarse %%%%

function [rhs, u0] = step_fine_to_coarse(Nx_pts, Nt_pts, L_h, R, rhs, u0, max_iter_sm, smoother, level, P)
%%%%%%%% loop to do pre smoothening

    %%% STEP 1: smooth L_h * u_approx = f by applying 
    % iterative scheme on fine grid
    if( level == 1)
     %fprintf('norm res before pre-smoothing : %d\n', norm(L_h*u0 - rhs));
    end
    
    if exist('P','var')
        u0 = smoothing(Nx_pts, Nt_pts, L_h, rhs, u0, max_iter_sm, smoother, P);
    else
        u0 = smoothing(Nx_pts, Nt_pts, L_h, rhs, u0, max_iter_sm, smoother);
    end

    if( level == 1)
     %res = L_h*u0 - rhs;
     %fprintf('norm res after pre-smoothing : %d\n', norm(res));
     %fprintf('norm res sigma before pre-smoothing : %d\n', norm(res(1:tot_pts)));
     %fprintf('norm res u pre-smoothing : %d\n', norm(res(tot_pts+1:end)));
    end
    
    %%%%% STEP 2: compute residual, then restrict residual
    % to coarse grid
    r = rhs - L_h*u0;
    
    % now restrict r{j}, and hence define new RHS
    rhs = R*r;    
end