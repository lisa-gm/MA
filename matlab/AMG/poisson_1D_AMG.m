%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% SOLVING A 1D POISSON EQUATION WITH AMG %%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% discretisation 
clear all;

N = 50; 
h = 1/(N-1);
x = linspace(0,1, N);

%%% here we know the real solution: 
u_real = 1/(4*pi^2) * sin(2*pi*x);

% and the RHS f
f = transpose(sin(2*pi*x));

dim = 1;
bdy_cond = 'Dirichlet';
L_h = set_up_L_h_FEM(N);
M_h = set_up_M_h_FEM(N);

% system we would like to solve L_h * u = f for u
u_init = zeros(N,1);
rhs = M_h*f;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AMG parameters
eps = 0.25;
max_iter = 100;
smoother = 'GaussSeidel';

tic;
u = two_level_AMG(L_h, rhs, u_init, smoother, eps, max_iter);
time = toc;

figure
plot(x, u_real, x, u);
legend('u real', 'u AMG');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot COARSE & FINE MESH

load('CoarseGridPts.mat');
C(C == 0) = NaN;
C(C == 1) = 0;
temp = zeros(length(x));

figure;
title('AMG Fine and Coarse Grid');
plot(x, temp, '.', x, C, '*');
legend('Fine Grid Points', 'Coarse Grid Points');
