%% setting up the RHS for DG in time and CG in space
% using P1 elements for both

function rhs = set_up_DG_rhs(dt, dx, Nt, Nx, f, M_t, M_h)
%% symbolic variables
syms t x

%% basis function on reference element q = 1 
% computing the integrals
% 
% l1(t) = 1-t/dt;
% l2(t) = t/dt;
% basis_t = {l1,l2};
% 
% k1(x) = 1-x/dx;
% k2(x) = x/dx;
% basis_x = {k1, k2};

%f_int = zeros(length(f),1);
% how to compute the integrals

% for k = 1:2
%     for l = 1:2
%         
%         bk = basis{k};
%         bl = basis{l};
%         
%         K_t(k,l) = - int(diff(bk,t)*bl,0,dt) + bk(dt)*bl(dt);              
%         
%         M_t(k,l) = int(bk*bl,0,dt);        
%         N_t(k,l) = bk(0)*bl(dt);       
%     end
% end

% for each time slab: get f_s, then putting them all together
% starting from f(0,x) ?!
f_matrix = reshape(f, [Nx, Nt]);
f_int_matrix = zeros(2*Nx, Nt-1);

for ts=1:Nt-1
    f_int_matrix(:, ts) = kron(M_t, M_h) * kron(f_matrix(:,ts), [1; 1]);
end

% reshape f into matrix with a column for each timestep
% so of size Nx x Nt
% duplicate columns
%f_dg = repmat(f_matrix, 2, 1);
% reshape f back into vector
%f_dg = reshape(f_dg, [Nx*Nt,1]);
% and then take out the first and last Nx entries...?!
rhs = reshape(f_int_matrix, [2*Nx*(Nt-1),1]);
end