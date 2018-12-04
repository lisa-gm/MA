function u = GaussSeidelSolve_LS(A, f, u0, max_iter)
tot_pts = length(f)/2;
% additional values, convergence threshold eps,
% damping factor
eps = 10^(-12);
w = 1;

for iter = 1:max_iter
    for i=1:tot_pts
        c = A([i,i+tot_pts],[i,i+tot_pts]) \ (f([i, i+tot_pts]) - A([i,i+tot_pts], :)*u0);
        u0([i, i+tot_pts]) = u0([i, i+tot_pts]) + c;
    end
    
    if(norm(f-A*u0) < eps)
       break;
    end
    
end


% get real solution
% u_real = A \ f;

% get iteration matrix
% D = diag(diag(A));
% E = - tril(A, -1);
% DL = (D/w-E);
% 
% % first iterate
% % use additional damping factor, can be just 1
% u_k = u0 + DL \ (f - (A * u0));
% 
% for k=2:max_iter
%     if(norm(f-A*u0) < eps)
%         break;
%     end
% 
%     u_k = u0 + DL\(f - (A * u0));
%     % new u becomes old u and initialise u_k to zero for summation
%     u0 = u_k;    
% end

u = u0;
end