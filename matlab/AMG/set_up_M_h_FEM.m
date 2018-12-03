%%%%%% set up M from heat eqn, that is matrix representing the integrals
% of phi_i * phi_j (test functions)

function M_h = set_up_M_h_FEM(N, bdy_type)
% distance between spacial coordinates
h = 1/(N-1);

% for now only do 1D
% matrix has a factor of h/6 in front and then 4's on diag
% 1 on both off diagonals
pts = N;
%pts = N-1;
temp = ones(pts,1);
M_h = h/6 * spdiags([temp, 4*temp, temp], -1:1, pts, pts);

if(strcmp(bdy_type, 'Dirichlet'))
    M_h(1,:) = [1, zeros(1,N-1)];
    M_h(end,:) = [zeros(1,N-1), 1];
end

% also eliminate entry in L_h(2,1) and L_h(N, N+1) to keep symmetry
% hence we have to update f(2) and f(N) later!
%M_h(2,1) = 0;
%M_h(N,N+1) = 0;

end