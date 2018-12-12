%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% smoother %%%%%%%%%%%%%%%%%%%%%

function u = smoothing(A, f, u, max_iter_sm, smoother)
    if(strcmp(smoother, 'Jacobi'))
        [u, ~] = JacobiSolve(A, f, u, max_iter_sm);
   
    elseif(strcmp(smoother, 'GaussSeidel'))
        u = GaussSeidelSolve(A, f, u, max_iter_sm);
            
    elseif(strcmp(smoother, 'Jacobi_LS'))
        u = JacobiSolve_LS(A, f, u, max_iter_sm);
        
    elseif(strcmp(smoother, 'Jacobi_LS_extended'))
        u = Jacobi_LS_extended(Nx_pts, Nt_pts, A, f, u0, max_iter, omega, block_size);     
                
    elseif(strcmp(smoother, 'GaussSeidel_LS'))
        u = GaussSeidelSolve_LS(A, f, u, max_iter_sm);    
    
    elseif(strcmp(smoother, 'ConjugateGradient'))
        [u, ~] = CGSolve(A, f, u, max_iter_sm);
        
    elseif(strcmp(smoother, 'GMRES'))
       u = gmres(A, f, max_iter_sm, 10^(-6), 1, eye(size(A)), eye(size(A)), u);
       
    else
        fprintf('invalid smoother!');
    end
    
end