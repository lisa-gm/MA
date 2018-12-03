%%%%%%%%%%%%%%%%%%%% ASSEMBLY SPACE TIME MATRIX %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% AND SUBMATRICES %%%%%%%%%%%%%%%%%%%%%%%

% function space : sum_l sum_k v_lk * phi(t) * psi(x) for 0<= k,l <= 1

% (K_t)_ij = \phi_i(t_m+1^-) * \phi_j(t_m+1^-) - int_over_t(\phi_i(t)*\phi'(t))dt
% (M_t)_ij = int_over_t(\phi_i(t)*\phi(t))dt
% (N_t)_ij = \phi_i(t_m^-) * \phi_j(t_m^+)
function [comp_mat, L_h, K_t, M_h, M_t, N_t] = set_up_M_st(ht, Nx, Nt, c, bdy_type)


[K_t, M_t, N_t] = set_up_DG_stencils(ht);

% Build the big matrices for entire space within one time slab
% (K_h)_ij = int_over_x( grad(psi_i)^T * grad(psi_j)) dx
% SCALE K_h BY CONSTANT C 
% double check that 1's get back in there where needed
L_h = c*set_up_L_h_FEM(Nx, bdy_type);

% (M_h)_ij = int_over_x( psi_i * psi_j) dx
M_h = set_up_M_h_FEM(Nx, bdy_type); 
 

%% build linear system
% this is for all elements of one time interval

A = kron(K_t, M_h) + kron(M_t, L_h);
B = -kron(N_t, M_h);

size_blk = size(A,1);

%% putting in identity rows
% %% Dirichlet BC
% doesnt seem to have an impact now but what if f \neq 0? 

if(strcmp(bdy_type,'Dirichlet'))
bdy_sides = [1 Nx Nx+1 2*Nx]; 

A(bdy_sides,:) = 0;
B(bdy_sides,:) = 0;

    for i = bdy_sides
        A(i,i) = 1;
    end
end

%% assemble

A_list = {};
B_list = {};
    
for i=1:Nt-1
    A_list{i} = A;
    B_list{i} = B;
end

A_list{1}(1:Nx, 1:Nx) = eye(Nx, Nx);
A_list{1}(1:Nx, Nx+1:end) = 0;

A_block_diag = sparse(blkdiag(A_list{1:end}));
B_block_low_diag = sparse(horzcat([zeros(size_blk, (Nt-2)*size_blk); blkdiag(B_list{1:end-1})], zeros((Nt-1)*size_blk, size_blk)));
    
comp_mat = A_block_diag + B_block_low_diag;

end
