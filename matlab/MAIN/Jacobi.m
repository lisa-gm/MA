%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% JACOBI SOLVE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function u = Jacobi(A, f, u_init, max_iter)
u = u_init;
w = 1;
eps = 10^-9;

invJ = diag(1./diag(A));

for i=1:max_iter
    u = u + w*invJ*(f - A*u);
    
    if(norm(f - A*u) < eps)
        % fprintf('converged after %d iterations\n', i);
        break;
    end
end

end