%% Trying to find out about the eigenvalues 

% suppose we choose square in time & space 
% of size 5 x 5, eg. time steps 3 - 7 and space elements 3 - 7
% which matrix entries does this correspond to? 
% can this be independent of rhs? 

% always the diagonal entries, 
% take entries: (3*10+3,3*10+3), (3*10+4,3*10+4), ... 

sp_ind = [3,4];
t_ind = [5,6];

indices = zeros(length(sp_ind)*length(t_ind),1);



%% setting up the problem parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nx = 10;
hx = 1/(Nx-1);
x = 0:hx:1;

T = 1;
Nt = 10;
ht = T/(Nt-1);
t = 0:ht:T;

k=1;
for i=1:length(sp_ind)
    for j=1:length(t_ind)
        indices(k) = Nx*(t_ind(j)-1)+sp_ind(i);
        k = k+1;
    end
end

% u0 = sin(pi*x)';
% if I use constant value
% have Pietros function for analytical solution
u0 = 5*ones(Nx,1);

f = zeros(Nx*Nt, 1);

bdy_left = zeros(Nt, 1);
bdy_right = zeros(Nt, 1);

u0 = [bdy_left(1); u0(2:end-1); bdy_right(1)];

% how do we handle the boundary points?
% each first and last row of each block e_i

%% the matrices and rhs

% function space : sum_l sum_k v_lk * phi(t) * psi(x) for 0<= k,l <= 1

% (K_t)_ij = \phi_i(t_m+1^-) * \phi_j(t_m+1^-) - int_over_t(\phi_i(t)*\phi'(t))dt
% (M_t)_ij = int_over_t(\phi_i(t)*\phi(t))dt
% (N_t)_ij = \phi_i(t_m^-) * \phi_j(t_m^+)

[K_t, M_t, N_t] = set_up_DG_stencils(ht);

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
B = -kron(N_t, M_h);

size_blk = size(A,1);

%% putting in identity rows
% %% Dirichlet BC
% doesnt seem to have an impact now but what if f \neq 0? 
bdy_sides = [1 Nx Nx+1 2*Nx]; 

A(bdy_sides,:) = 0;
B(bdy_sides,:) = 0;

for i = bdy_sides
    A(i,i) = 1;
end

%% assemble

A_list = {};
B_list = {};
    
for i=1:Nt-1
    A_list{i} = A;
    B_list{i} = B;
end
    
A_block_diag = sparse(blkdiag(A_list{1:end}));
B_block_low_diag = sparse(horzcat([zeros(size_blk, (Nt-2)*size_blk); blkdiag(B_list{1:end-1})], zeros((Nt-1)*size_blk, size_blk)));
    
comp_mat = A_block_diag + B_block_low_diag;

eig_val = eig(full(comp_mat(indices, indices)))

%% computing f properly
% need int_over_all(f(x,t)*phi(t)*psi(x))dxdt
rhs = set_up_DG_rhs(ht, hx, Nt, Nx, f, M_t, M_h);

%%%%%%% OPTION 3 %%%%%%%

comp_mat_3 = comp_mat;

rhs_3 = rhs;
% is this not going to work anymore if f \neq 0 ?

u0_f = kron([1,0],u0');  
rhs_3(1:2*Nx) = kron([1 0 ; 0 0],M_h) * u0_f';

%% solve
% comp_mat * u = rhs

u_3 = comp_mat_3 \ rhs_3;

%% rearrange

u_mat_3 = zeros(Nx, Nt);

u_mat_3(:,1) = u0;

% note that we are only taking one set of solutions for 
% the inner nodes, always x_i- 
% WHY better than x_i+ ?? 

for ts=1:Nt-2
    u_mat_3(:, ts+1) = u_3((2*ts-1)*Nx+1:2*ts*Nx);
end

%% plot

% figure
% mesh(x, t, u_mat_1');
% title('approximated solution');


%% computing the analytical solution
% for u0 constant value

u0 = 5;

% compute analitical solution
u_an = zeros(size(u_mat_1));
%u_true = zeros(1,length(u));
x = linspace(0,1,Nx);
t = linspace(0,T,Nt);

for i = 1:length(x)
   for j = 1:length(t)       
      u_an(i,j) = heat_1D_an( u0,x(i),t(j),Nt);        
   end 
end

u_an(1:Nx,1) = [ 0; u0*ones(Nx-2,1);0]; % initial condition

% figure
% mesh(x,t, u_an');
% title('analytical solution')

%% comparison with analytical solution

error_3 = zeros(1, Nt);


for i=1:Nt

    error_3(i) = norm(u_mat_3(:,i) - u_an(:,i));

end

% figure
% plot(0:Nt-1, error_3);
% legend('opt 3');
% title('Error compared to analytical solution each timestep');