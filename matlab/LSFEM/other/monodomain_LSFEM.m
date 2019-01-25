%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% NEWTON SOLVER FOR MONODOMAIN %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;

%% set up general parameters

max_iter_n = 1;
eps = 10^(-6);

%% problem parameters 

%%%%%%%%%%%%%%%%%%%%%%%% GRID %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nx = 20;
S = 1;
hx = S/(Nx-1);
x = linspace(0, S, Nx);

Nt = 80;
T = 8;
ht = T/(Nt-1);
t = linspace(0, T, Nt);

N_tot = Nx*Nt;

%%%%%%%%%%%%%%  FORCING TERM F %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f = @(u) u.*(1-u).*(u-0.1);
deriv_f = @(u) -3*u.*u + 2.2*u - 0.1;

%%%%%%%%%%%%% INITIAL AND BOUNDARY CONDITIONS %%%%%%%%%%%%%%%%%

% initial conditions at t=0
u0 = zeros(Nx,1);
for i=1:Nx
    u0(i) = max(0, 1-2*x(i));
end
% u0 = 5*ones(Nx,1);

%u0 = transpose(-2.5*x.^2 + 2.5*x+1);
%u0 = transpose(1/(pi^2)*sin(pi*x)+1);
%u0 = zeros(Nx, 1);

% BOUNDARY VALUE IN TIME 

%bdy_left = ones(Nt, 1);
%bdy_right = ones(Nt, 1);
bdy_left = zeros(Nt, 1);
bdy_right = zeros(Nt, 1);
%bdy_left = linspace(0,1, Nt);
%bdy_right = linspace(0,1,Nt);

%%%%%%%%%%%%%%% INITIAL GUESS FOR u(x,t) %%%%%%%%%%%%%%%%%

u = 0.1*ones(N_tot, 1);
u(1:Nx) = u0;
u(1:Nx:end) = bdy_left;
u(Nx:Nx:end) = bdy_right;

%% constructing the operators

[G_h, K_h, M_h] = set_up_FEM_1D_mat(Nx);
[G_t, K_t, M_t] = set_up_FEM_1D_mat(Nt);

%%%%%%%%%%%%%%%%%%% LINEAR PART %%%%%%%%%%%%%%%%%%%%%%%%%%
% ASSEMBLE A

A_ss = T*S*kron(M_t, M_h) + T/S*kron(M_t, K_h);
A_su =  -T*kron(M_t, G_h) - kron(G_t, G_h');
A_us = - T*kron(M_t, G_h') - kron(G_t', G_h);
A_uu = + T/S*kron(M_t, K_h) + S/T*kron(K_t, M_h);

%%%%%%%%%%%%%%%%%% NON LINEAR PART %%%%%%%%%%%%%%%%%%%%%%%
% has to be done in loop

%% loop






%% plotting
