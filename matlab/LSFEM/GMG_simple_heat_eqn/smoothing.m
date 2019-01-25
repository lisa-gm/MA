%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% smoother %%%%%%%%%%%%%%%%%%%%%

function u = smoothing(Nx_pts, Nt_pts, A, f, u, max_iter_sm, smoother, P)

if exist('P','var')
     omega = 0.64;
     [u, ~] = apply_smoother(A, P, f, u, omega, max_iter_sm);
    
 
else
    if(strcmp(smoother, 'Jacobi'))
        [u, ~] = JacobiSolve(A, f, u, max_iter_sm);
        
    elseif(strcmp(smoother, 'line_smoother'))
        u = line_smoother(Nx_pts, Nt_pts, A, f, u, max_iter_sm, 10^(-12));             
   
    elseif(strcmp(smoother, 'GaussSeidel'))
        u = GaussSeidelSolve(A, f, u, max_iter_sm);
        
    elseif(strcmp(smoother, 'two_block_GaussSeidel'))
        u = two_block_GaussSeidel(A, f, u, max_iter_sm, 10^(-12));   
            
    elseif(strcmp(smoother, 'Jacobi_LS'))
        u = JacobiSolve_LS(A, f, u, max_iter_sm);
        
    elseif(strcmp(smoother, 'Jacobi_LS_extended'))
        omega = 0.63;
        block_size = 2;
        u = JacobiSolve_LS_extended(Nx_pts, Nt_pts, A, f, u, max_iter_sm, omega, block_size);     
    
    elseif(strcmp(smoother, 'Jacobi_LS_extended_SP_T'))
        omega = 0.64;
        block_size_s = 2;
        block_size_t = 2;
        u = JacobiSolve_LS_extended_SP_T(Nx_pts, Nt_pts, A, f, u, max_iter_sm, omega, block_size_s, block_size_t);                 

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
    
end
