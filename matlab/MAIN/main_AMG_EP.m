%%%%%%%%%%%%%%% ALGEBRAIC MULTIGRID %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2D + time 
% order of the set up: first walk in x-direction, then 
% in y-direction, then time.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% how do I incorporate RIGHT HAND SIDE ?!? i
% just use f depending on previous values of sol V ?!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% set up of fine grid %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% space-time version, Ganda & Neumueller approach ?!
% use 5 pt stencil for space and implicit Euler 
% or is an upwind scheme better ?! 

% two space dimensions, for now just square [0,1]^2
Nx = 10;
x = linspace(0,1, Nx);
hx = 1/(Nx-1);

N_sp = Nx*Nx;

% time steps
Nt = 10;
T = 1;
t = linspace(0,T, Nt);
ht = T/(Nt-1);

N_tot = N_sp * Nt;

% initial conditions first as matrix
V0 = zeros(Nx, Nx);
mid = floor(Nx/2);
V0(mid-1:mid+1, mid-1:mid+1) = 1;

% boundary conditions space
% assume for now constant over time
temp = zeros(Nx, 1);

% there is an overlap, need to make sure, they are the same
bdy_north = temp;
bdy_south = temp;
bdy_west = temp;
bdy_east = temp;

% set up RHS
% keep recalling them

% FEM space operator, stiffness matrix, 1 timestep
Lh = set_up_Lh_FEM_2D(Nx);

% FEM space operator, mass matrix, 1 timestep
Mh = set_up_Mh_FEM_2D(Nx);

% DG in time
[K_t, M_t, N_t] = set_up_DG_stencils(ht);

% add in boundary conditions into overall solution vector
V_mat = zeros(Nx, Nx, Nt);

% TODO: get rid of loop, more elegant way ...
for i=1:Nx
    V_mat(1, :, i) = bdy_south;
    V_mat(Nx, :, i) = bdy_north;
    V_mat(:, 1, i) = bdy_west;
    V_mat(:, Nx, i) = bdy_east;
end

%% set up overall space time matrix %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A = kron(K_t, M_h) + kron(M_t,K_h);
B = - kron(N_t, M_h);


tot_mat = 



%% monodomain model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dV/dt - delta(V) = I_ext + I_ion
% where I_ext =  
% I_ion = V * (a - V) * (V - 1)

%% WHAT DO I NEED: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% interpolation + restriction operator
% coarse level A

% how to compute the weights?


%% want to analyse which are the taken coarse grid points 
% and the respective weights


