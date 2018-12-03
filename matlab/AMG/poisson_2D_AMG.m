%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% SOLVING A 2D POISSON EQUATION WITH AMG %%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% discretisation 
clear all;

tic; 
Nx = 10; 
hx = 1/(Nx -1);
x = linspace(0,1, Nx);

Ny = 12; 
hy = 1/(Ny -1);
y = linspace(0,1, Ny);

%%% here we know the real solution: 
%u_real = 1/(4*pi^2) * sin(2*pi*x);

% and the RHS f, what would be a good RHS ?!
f = zeros(Nx, Ny);

for j=1:Ny
    for i=1:Nx
    f(i,j) = sin(2*pi*x(i))*cos(2*pi*y(j));
    end
end

f_vec = reshape(f, [Nx*Ny, 1]);

dim = 2;
bdy_cond = 'Dirichlet';
L_h = set_up_L_h_FEM_2D(Nx, Ny, bdy_cond);
M_h = set_up_M_h_FEM_2D(Nx, Ny, bdy_cond);

% system we would like to solve L_h * u = f for u
u_init = 1/(4*pi^2)*sin(2*pi*x);
rhs = M_h*f_vec;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% incorporate initial conditions
L_h(1:Nx, :) = horzcat(eye(Nx, Nx), zeros(Nx, (Ny-1)*Nx));
rhs(1:Nx) = u_init';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% direct solve

% u_large = L_h \ rhs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AMG parameters
eps = 0.25;
max_iter = 100;
smoother = 'GaussSeidel';
u_guess = [u_init'; 0.1*ones(Nx*(Ny-1),1)];

tic;
u_large = two_level_AMG(L_h, rhs, u_guess, smoother, eps, max_iter);
AMG_time = toc;

fprintf('time AMG: %d sec\n', time);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reshape u  &  PLOTTING
 
u_mat = zeros(Nx, Ny);

for s=1:Ny
    u_mat(:, s) = u_large((s-1)*Nx+1:s*Nx);
end

figure
mesh(x, y, u_mat')
title('2D Poisson Eqn');
xlabel('x axis');
ylabel('y axis');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot COARSE & FINE MESH

load('CoarseGridPts.mat');
C(C == 0) = NaN;
C(C == 1) = 0;

F_mat = zeros(Nx, Ny);
C_mat = ones(Nx, Ny);
for s=1:Ny
    C_mat(:, s) = C((s-1)*Nx+1:s*Nx);
end

[X, Y] = meshgrid(x,y);

figure;
plot3(X, Y, F_mat', 'b.', X, Y, C_mat', 'r*');
legend('fine level gridpts', 'coarse level grid pts');
title('Fine and Coarse Level Grid Points');
xlabel('x axis');
ylabel('y axis');

overall_time = toc;
fprintf('overall time: %d sec \n', overall_time);
