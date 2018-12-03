%% setting up the RHS for DG in time and CG in space
% using P1 elements for both

function rhs = set_up_DG_rhs(dt, dx, Nt, Nx, f, M_t, M_h)

% for each time slab: get f_s, then putting them all together
% starting from f(0,x) ?!
f_matrix = reshape(f, [Nx, Nt]);
f_int_matrix = zeros(2*Nx, Nt-1);

for ts=1:Nt-1
    f_int_matrix(:, ts) = kron(M_t, M_h) * kron(f_matrix(:,ts), [1; 1]);
end

% get into vector format
rhs = reshape(f_int_matrix, [2*Nx*(Nt-1),1]);
end