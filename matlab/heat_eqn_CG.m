%%% Heat Equation 
% solving with CG to compare with DG

% also use Ganda + Neumueller approach to discretise in space time
% one term less from weak form

% A = K_tilde_t kron M_h + M_t kron K_h
% B = 0

%% setting up the problem parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nx = 10;
hx = 1/(Nx-1);
x = 0:hx:1;

T = 1;
Nt = 10;
ht = T/(Nt-1);
t = 0:ht:T;

% u0 = sin(pi*x)';
% if I use constant value
% have Pietros function for analytical solution
u0 = 5*ones(Nx,1);

f = zeros(Nx*Nt, 1);

bdy_left = zeros(Nt, 1);
bdy_right = zeros(Nt, 1);

u0 = [bdy_left(1); u0(2:end-1); bdy_right(1)];

%% set up stencils

[K_t, M_t] = set_up_CG_stencils_time(ht);

% Build the big matrices for entire space within one time slab
% (K_h)_ij = int_over_x( grad(psi_i)^T * grad(psi_j)) dx
K_h = set_up_L_h_FEM(Nx);
% (M_h)_ij = int_over_x( psi_i * psi_j) dx
M_h = set_up_M_h_FEM(Nx); 

% at the end we want ones on the at the beginning and end of each block
% therefore we set M_h to zero 
M_h(1,1) = 0;
M_h(end,end) = 0;

%% build linear system
% this is for all elements of one time interval

A = kron(K_t, M_h) + kron(M_t, K_h);

size_blk = size(A,1);

%% putting in identity rows
% %% Dirichlet BC
% doesnt seem to have an impact now but what if f \neq 0? 
bdy_sides = [1 Nx Nx+1 2*Nx]; 

A(bdy_sides,:) = 0;

for i = bdy_sides
    A(i,i) = 1;
end

%% assemble

A_list = {};
    
for i=1:Nt-1
    A_list{i} = A;
end
    
comp_mat = sparse(blkdiag(A_list{1:end}));

% need int_over_all(f(x,t)*phi(t)*psi(x))dxdt

rhs = set_up_DG_rhs(ht, hx, Nt, Nx, f, M_t, M_h);
   