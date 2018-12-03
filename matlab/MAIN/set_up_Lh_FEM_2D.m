%% regular 5-point stencil Continuous Galerkin %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% on Lh_ii = 4, on Lh_ii-1 = Lh_ii+1 = Lh_i-1i = Lhi+1i = -1
% otherwise zero
% add in boundary conditions later

%%%%% USING A TRIANGULAR ELEMENTS %%%%%%%%%%%%%%%%%

function Lh = set_up_Lh_FEM_2D(Nx)

hx = 1/(Nx-1);
N_tot = Nx*Nx;

%%% setting up the diagonal block matrices %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

e_small = ones(Nx, 1);
block = spdiags([-1*e_small, 4*e_small, -1*e_small], -1:1, Nx, Nx);

block_list = {};

for i=1:Nx
    block_list{i} = block;
end

Lh_diag = blkdiag(block_list{1:end})

%%%%%%% setting up the off diagonal blocks %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

e_large = ones(N_tot, 1);
Lh_sides = spdiags([-1*e_large, -1*e_large], [-Nx, Nx], N_tot, N_tot);

%%%%% adding them together + scaling %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Lh = 1/hx^2.*(Lh_sides + Lh_diag);

%%%% putting in zeros and ones for boundary %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Lh(1:Nx:end,:) = zeros(Nx, N_tot);
Lh(Nx:Nx:end,:) = zeros(Nx, N_tot);

Lh(1:Nx:end, 1:Nx:end) = eye(Nx);
Lh(Nx:Nx:end, Nx:Nx:end) = eye(Nx);

end
