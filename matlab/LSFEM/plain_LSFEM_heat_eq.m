%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% SETTING UP LEAST SQUARES SPACE TIME HEAT EQ %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% u_t - delta(u) = f(x,t) %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%% u \in P1P1 AND \sigma in P1P1 %%%%%%%%%%%%%%%%%%%%%
clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% GENERAL PARAMETERS 

Nx = 20;
S = 1;
hx = S/(Nx-1);
x = linspace(0, S, Nx);

Nt = 20;
T = 1;
ht = T/(Nt-1);
t = linspace(0, T, Nt);

%%%%%%%%%%%% initial conditions at t = 0 %%%%%%%%%%%%%%
%u0 = transpose(-2.5*x.^2 + 2.5*x+1);
%u0 = transpose(1/(pi^2)*sin(pi*x)+1);
%u0 = ones(Nx,1);
u0 = linspace(0, S, Nx);

%%%%%%%%%%%% setting up rhs f %%%%%%%%%%%%%%%%%%%%%%%%%
f = zeros(Nx, Nt);
f_int_sigma = zeros(size(f));
f_int_u = zeros(size(f));

for ts=1:Nt
    %f(:,ts) = 5;
    f(:,ts) = sin(pi*x)*cos(pi*t(ts)); %+ cos(2*pi*t(ts));

end

f_vec = reshape(f, [(Nx*Nt),1]);

%%%%%%%%% BOUNDARY CONDITIONS in TIME %%%%%%%%%%%%%%%%%%%%%%%%%%%%

%bdy_left = ones(Nt, 1);
bdy_right = ones(Nt, 1);
bdy_left = zeros(Nt, 1);
%bdy_right = zeros(Nt, 1);
%bdy_left = linspace(0,1, Nt);
%bdy_right = linspace(0,1, Nt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% GET STIFFNESS, GRADIENT, MASS Matrices in time + space

[G_h, K_h, M_h] = set_up_FEM_1D_mat(Nx);
[G_t, K_t, M_t] = set_up_FEM_1D_mat(Nt);

% ASSEMBLE A

A_ss = kron(M_t, M_h) + kron(M_t, K_h);
A_su =  kron(M_t, G_h) + kron(G_t', G_h);
A_us = + kron(M_t, G_h') + kron(G_t, G_h');
A_uu = + kron(M_t, K_h) + kron(K_t, M_h);

% check for symmetry: norm(full(A_us - A_su'))

% ENFORCING BOUNDARY CONDITIONS

A_uu(1:Nx:end, :) = 0;
A_uu(Nx:Nx:end, :) = 0;
A_uu(1:Nx, :) = 0;

A_uu(1:Nx, 1:Nx) = eye(Nx);
A_uu(1:Nx:end, 1:Nx:end) = eye(Nt);
A_uu(Nx:Nx:end, Nx:Nx:end) = eye(Nt);

A_us(1:Nx, :) = 0;
A_us(1:Nx:end, :) = 0;
A_us(Nx:Nx:end, :) = 0;

%%%%%%%%%%%% ASSEMBLING A %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A = [A_ss, A_su; A_us, A_uu];

size_block = size(A_ss, 1);
size_A = size(A,1);

% ASSEMBLE RHS

f_s = kron(M_t, G_h') * f_vec;
f_u = kron(G_t', M_h) * f_vec;

f_u(1:Nx) = u0;
f_u(1:Nx:end) = bdy_left;
f_u(Nx:Nx:end) = bdy_right;

F = [f_s; f_u];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% SOLVE

sol = A \ F;
u = sol(size_block+1:end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% PLOT

% reshape
u_mat = zeros(Nx, Nt);

for ts=1:Nt
    u_mat(:, ts) = u((ts-1)*Nx+1:ts*Nx);
end

figure
mesh(x, t(1:end), transpose(u_mat(:,1:end)));
xlabel('space');
ylabel('time');
zlabel('u');
zlim([min(u), max(u)]);
title('heat equation LSFEM');

