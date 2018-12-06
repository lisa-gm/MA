function u = GaussSeidelSolve_LS(A, f, u0, max_iter)
tot_pts = length(f)/2;
% additional values, convergence threshold eps,
% damping factor
eps = 10^(-12);

for iter = 1:max_iter
    for i=1:tot_pts
        %tot_pts
        %size(A)
        %A([i,i+tot_pts],[i,i+tot_pts])
        %f([i, i+tot_pts])
        c = A([i,i+tot_pts],[i,i+tot_pts]) \ (f([i, i+tot_pts]) - A([i,i+tot_pts], :)*u0);
        u0([i, i+tot_pts]) = u0([i, i+tot_pts]) + c;
    end
    
    if(norm(f-A*u0) < eps)
       break;
    end
    
    %fprintf('norm res : %d\n', norm(f-A*u0));
    
end

u = u0;
end