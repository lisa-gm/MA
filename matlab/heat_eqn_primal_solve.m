%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% SETTING UP LEAST SQUARES SPACE TIME HEAT EQ %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%% u \in P1P1 AND \sigma in P0P1 %%%%%%%%%%%%%%%%%%%%%
%clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% GENERAL PARAMETERS 

Nx = 40;
S = 1;
hx = S/(Nx-1);
x = linspace(0, S, Nx);

Nt = 500;
T = 10;
ht = T/(Nt-1);
t = linspace(0, T, Nt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% FORCING TERM / RHS 

f = zeros(Nx, Nt);
%f = ones(Nx, Nt);

for ts=1:Nt
    %f(:,ts) = 5;
    f(:,ts) = sin(pi*x)*cos(pi*t(ts)); %+ cos(2*pi*t(ts));
   
end

f_vec = reshape(f, [(Nx*Nt),1]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% INITIAL AND BOUNDARY CONDITIONS

%u0 = transpose(-2.5*x.^2 + 2.5*x+1);
u0 = transpose(1/(pi^2)*sin(pi*x)+1);
%u0 = zeros(Nx, 1);
%u0 = x';

bdy_left = ones(Nt, 1);
bdy_right = ones(Nt, 1);

%bdy_left = zeros(Nt, 1);
%bdy_right = zeros(Nt, 1);
%bdy_left = linspace(0,1, Nt);
%bdy_right = linspace(0,1,Nt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% GET STIFFNESS, GRADIENT, MASS Matrices in time + space

[G_h, K_h, M_h] = set_up_FEM_1D_mat(Nx);
[G_t, K_t, M_t] = set_up_FEM_1D_mat(Nt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%  PRIMAL SPACE TIME FORMULATION FEM 

% scaling to compare with other programs, cancels in this formulation

scale = 1/(ht*hx);
%scale = 1;

A_primal = full(S*scale*kron(G_t, M_h) + T/S*scale*kron(M_t, K_h)); %
%A_primal = full(S*scale*kron(G_t, M_h));

A_primal(1:Nx,:) = 0;
A_primal(1:Nx:end, :) = 0;
A_primal(Nx:Nx:end, :) = 0;

A_primal(1:Nx, 1:Nx) = eye(Nx);
A_primal(1:Nx:end, 1:Nx:end) = eye(Nt);
A_primal(Nx:Nx:end, Nx:Nx:end) = eye(Nt);

f_primal = T*S*scale*kron(M_t, M_h)*f_vec;

f_primal(1:Nx) = u0;
f_primal(1:Nx:end) = bdy_left;
f_primal(Nx:Nx:end) = bdy_right;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% SOLVE

u_primal = A_primal\f_primal;

% reshape solution
u_mat_prim = zeros(Nx, Nt);

for ts=1:Nt
    u_mat_prim(:, ts) = u_primal((ts-1)*Nx+1:ts*Nx);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% PLOT

figure
mesh(x, t(1:end), transpose(u_mat_prim(:,1:end)));
xlabel('space');
ylabel('time');
zlabel('u');
title('heat equation from primal weak formulation');

