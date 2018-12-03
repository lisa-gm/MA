%% used construct the discretised weak form of any vector f
% with Nx space components per time slab and Nt time slabs
% DOES NOT CHANGE SIZE
% ie if already in DG w/ double entries stays like that 
% if not then also stays like that

function rhs = set_up_mass_mat_rhs(dt, dx, slices, Nx, f, M_t, M_h)

f_matrix = reshape(f, [Nx, slices]);
f_int_matrix = zeros(size(f_matrix));

for ts=1:slices
    f_int_matrix(:, ts) = kron(M_t, M_h) * f_matrix(:,ts);
end

% get into vector format
rhs = reshape(f_int_matrix, [Nx*slices,1]);
end