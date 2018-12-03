%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% SETTING UP LEAST SQUARES POISSON EQ %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% GENERAL PARAMETERS 

Nx = 5;
S = 1;
hx = S/(Nx-1);
x = linspace(0, S, Nx);

f_int = zeros(Nx,1);

for i=1:Nx-1
    tmp=1/hx*int(f, y, x(i), x(i+1));
    f_int(i) = f_int(i)-tmp;
    f_int(i+1)=f_int(i+1)+tmp;
end

%f_vec = 5*ones(Nx,1);
f_vec = -transpose(sin(pi*x));
%f_vec = transpose(2*x.^2);

%u_real = 2.5*x.^2 -2.5*x; % u_real' = 5*x - 2.5 --> at 1: 
%sigma_left_dir = 2.5;
%sigma_right_dir = -2.5;

% bdy conditions
bdy_left = 1;
bdy_right = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% GET STIFFNESS, GRADIENT, MASS Matrices in space

[G_h, K_h, M_h] = set_up_FEM_1D_mat(Nx);

% ASSEMBLE A
% took out the terms that depended on d/dt, doesnt exist here

A_ss = M_h + K_h;
A_su = - G_h; 
A_us = - G_h'; 
A_uu = K_h;

A = [A_ss, A_su; A_us, A_uu];

size_block = size(A_ss, 1);
size_A = size(A,1);

A(size_block+1, 1:end) = [zeros(1, size_block), 1, zeros(1, size_block-1)];
A(end, 1:end) = [zeros(1, size_A-1), 1];

%A(1, 1:end) = [1, zeros(1, 2*size_block-1)];
%A(size_block, 1:end) = [zeros(1, size_block-1), 1, zeros(1,size_block)];

% ASSEMBLE RHS

f_s = -G_h'*f_vec;
% trying to impose conditions on f_s
%f_s(1) = sigma_left_dir;
%f_s(end) = sigma_right_dir;

f_u = zeros(Nx,1);   % since dv/dt = 0 ....
f_u(1) = bdy_left; f_u(end) = bdy_right;

F = [f_s; f_u];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% SOLVE

sol = A \ F;
u = sol(size_block+1:end);

% real solution
K_h(1, 1:end) = [1, zeros(1, Nx-1)];
K_h(end, 1:end) = [zeros(1, Nx-1), 1];

rhs = M_h*f_vec;
rhs(1) = bdy_left; rhs(end) = bdy_right;

u_real = K_h \ rhs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% PLOT 

% ,x, 1/(4*pi^2)*sin(2*pi*x) , '1/(4*pi^2)*sin(2*pi*x)'

figure;
plot(x, u, x, u_real);
title('Poisson Eq solved w/ LSFEM');
legend('LSFEM sol', 'FEM backslash sol');
xlabel('space');
