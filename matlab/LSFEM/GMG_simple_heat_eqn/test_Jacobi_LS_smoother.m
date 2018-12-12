% ******************************************************************* %
% *********************** TEST JACOBI SOLVE ************************* %
% ******************************************************************* %


max_iter = 5;

Nx_pts = 2;
Nt_pts = 2;

f = rand(2*Nx_pts*Nt_pts,1);
A = rand(2*Nx_pts*Nt_pts);
u0 = zeros(size(f));

u_direct = A \ f;

[u, sol] = JacobiSolve_LS(A, f, u0, max_iter);

%[u_ex, sol_ex] = JacobiSolve_LS_extended(Nx_pts, Nt_pts, A, f, u0, max_iter, 0.6656, 2);