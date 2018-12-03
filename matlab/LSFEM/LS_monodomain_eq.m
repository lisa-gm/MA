%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% SETTING UP LEAST SQUARES SPACE TIME MONODOMAIN %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% GENERAL PARAMETERS 

max_iter = 5;
c = 10^(-3);
a = 10;
omega = 1;

Nx = 40;
S = 1;
hx = S/(Nx-1);
x = linspace(0, S, Nx);

Nt = 10;
T = 1;
ht = T/(Nt-1);
t = linspace(0, T, Nt);

% INITIAL CONDITIONS, BDY CONDITIONS and F
% initial conditions at t=0
u0 = zeros(Nx,1);
for i=1:Nx
    u0(i) = max(0, 1-2*x(i));
end

f = @(u) a*u.*(1-u).*(u-0.1);
deriv_f = @(u) -3*a*u.*u + 2.2*a*u - a*0.1;

%%%%%%%%% BOUNDARY CONDITIONS

% homogenous Neumann, so nothing except for u0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% GET STIFFNESS, GRADIENT, MASS Matrices in time + space %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[G_h, K_h, M_h] = set_up_FEM_1D_mat(Nx);
[G_t_u, K_t_u, M_t_u] = set_up_FEM_1D_mat(Nt);

% change basis functions of \sigma, \tau to
% \sigma = sum_j \sigma_j*\psi_j bec. no time derivative there
% hence M_t = 1/ht * eye(Nt);
[G_p0p1, M_p0p1, M_p0p0] = set_up_P0P1_FEM_mat(Nt);

% ASSEMBLE A

A_ss = kron(M_p0p0, M_h) + kron(M_p0p0, c^2*K_h);
A_su = - kron(M_p0p1, G_h) - kron(G_p0p1, c*G_h);
A_us = - kron(M_p0p1, G_h') - kron(G_p0p1', c*G_h');
A_uu = kron(M_t_u, K_h) + kron(K_t_u, M_h);


% ENFORCING BOUNDARY CONDITIONS
% need to do this on rhs to keep symmetry ... 
%A_su(end-Nx+1:end, :) = 0;
%A_ss(end-Nx+1:end, end-2*Nx+1:end-Nx) = -eye(Nx);
%A_ss(end-Nx+1:end, end-Nx+1:end) = eye(Nx);

%A_uu(1:Nx:end, :) = 0;
%A_uu(Nx:Nx:end, :) = 0;
A_uu(1:Nx, :) = 0;
%A_uu(end-Nx+1:end, :) = 0;

A_uu(1:Nx, 1:Nx) = eye(Nx);
%A_uu(end-Nx+1:end, end-Nx+1:end) = eye(Nx);
%A_uu(end-Nx+1:end, end-2*Nx+1:end-Nx) = -eye(Nx);
%A_uu(1:Nx:end, 1:Nx:end) = eye(Nt);
%A_uu(Nx:Nx:end, Nx:Nx:end) = eye(Nt);

A_us(1:Nx, :) = 0;
%A_us(1:Nx:end, :) = 0;
%A_us(Nx:Nx:end, :) = 0;
%A_us(end-Nx+1:end, :) = 0;

A = [A_ss, A_su; A_us, A_uu];

size_block = size(A_ss, 1);
size_A = size(A,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% NEWTON ITERATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N_tot = size(A, 1);
sigma = 0.0*ones(size_block, 1);
u = 0.2*ones(size_block, 1);
u(1:Nx) = u0;

sol = [sigma; u];
sol_old = 10*ones(size_A, 1);

M_p0_G_h = kron(M_p0p0, c*G_h);
G_tu_M_h = kron(G_t_u, M_h);

iter = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% LOOP

% what is a good stopping criterion
while(iter <= max_iter && norm(sol - sol_old) > eps)
    tic;
    
    f_s = M_p0_G_h * f(u);
    f_u = G_tu_M_h * f(u);
    f_u(1:Nx) = u0;

    F = [f_s; f_u];

    J = A*sol - F;

    % just define first entries to be zero, its where IC are
    J(size_block+1:size_block+Nx) = zeros(Nx, 1);
    
    % construct Hessian: 
    % H = A - diag(M_p0_G_h*deriv_f(u), G_tu_M_h*deriv_f(u))
    
    der_f = deriv_f(u);
    H = sparse(A - diag([zeros(size_block,1); G_tu_M_h*der_f]));

    u_mat = zeros(Nx, Nt);

    for ts=1:Nt
    u_mat(:, ts) = u((ts-1)*Nx+1:ts*Nx);
    end

    figure
    mesh(x, t(1:end), transpose(u_mat(:,1:end)));
    xlabel('space');
    ylabel('time');
    zlabel('u');
    title('heat equation LSFEM');
    
    % check if H diagonally dominant
%     arg_max = zeros(1, size(H,1));
%     for row=1:size(H,1)
%         [~, arg_max(row)] = max(abs(H(row,:)));
%     end
%     
%     
%     temp = arg_max - 1:size(H,1);
%     temp(temp == 0) = 1;
%     temp(temp ~= 0) = 0;
    
    
    sol_old = sol;
    
    % direct solve
    sol = sol_old - omega*(H \ J);
    u = sol(size_block+1:end);
  
    iter = iter+1;
    tt = toc;
%     fprintf('one complete Newton step incl AMG, time: %d sec, iter: %d \n', tt, iter-1);
end

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
title('heat equation LSFEM');

