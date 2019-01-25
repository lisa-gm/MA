%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% SETTING UP LEAST SQUARES SPACE TIME HEAT EQ %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%% u \in P1P1 AND \sigma in P0P1 %%%%%%%%%%%%%%%%%%%%%
%clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% GENERAL PARAMETERS 

Nx = 10;
S = 1;
hx = S/(Nx-1);
x = linspace(0, S, Nx);

Nt = 40;
T = 4;
ht = T/(Nt-1);
t = linspace(0, T, Nt);

u0 = transpose(-2.5*x.^2 + 2.5*x+1);
%u0 = transpose(1/(pi^2)*sin(pi*x)+1);
%u0 = zeros(Nx, 1);

% f
%f = zeros(Nx, Nt);
f = ones(Nx, Nt);

for ts=1:Nt
    f(:,ts) = 5;
    %f(:,ts) = sin(pi*x)*cos(pi*t(ts)); %+ cos(2*pi*t(ts));
   
end

f_vec = reshape(f, [(Nx*Nt),1]);

%%%%%%%%% BOUNDARY CONDITIONS

bdy_left = ones(Nt, 1);
bdy_right = ones(Nt, 1);
%bdy_left = zeros(Nt, 1);
%bdy_right = zeros(Nt, 1);
%bdy_left = linspace(0,1, Nt);
%bdy_right = linspace(0,1,Nt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% GET STIFFNESS, GRADIENT, MASS Matrices in time + space

[G_h, K_h, M_h] = set_up_FEM_1D_mat(Nx);
[G_t_u, K_t_u, M_t_u] = set_up_FEM_1D_mat(Nt);
[G_p0p1, M_p0p1, M_p0p0] = set_up_P0P1_FEM_mat(Nt);

% ASSEMBLE A

A_ss = kron(M_p0p0, M_h) + kron(M_p0p0, K_h);
A_su = - kron(M_p0p1, G_h) - kron(G_p0p1, G_h');
A_us = - kron(M_p0p1', G_h') - kron(G_p0p1', G_h);
A_uu = + kron(M_t_u, K_h) + kron(K_t_u, M_h);

norm(full(A_us - A_su'))
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

A = [A_ss, A_su; A_us, A_uu];

size_block = size(A_ss, 1);
size_A = size(A,1);

% ASSEMBLE RHS
% T0 DO: deal with f_s appropriately 
f_s = -kron(M_p0p0, G_h') * f_vec(Nx+1:end);
f_u = kron(G_t_u', M_h) * f_vec;

f_u(1:Nx) = u0;

f_u(1:Nx:end) = bdy_left;
f_u(Nx:Nx:end) = bdy_right;

F=  [f_s; f_u];

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%% SOLVE

sol = A \ F;
u = sol(size_block+1:end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% PLOT

% reshape
u_mat = zeros(Nx, Nt);


for ts=1:Nt
    % u_mat(:, ts) = u((ts-1)*Nx+1:ts*Nx);
    u_mat(:, ts) = u((ts-1)*Nx+1:ts*Nx);
end

figure
mesh(x, t(1:end), transpose(u_mat(:,1:end)));
xlabel('space');
ylabel('time');
zlabel('u');
title('heat equation LSFEM P0P1');


