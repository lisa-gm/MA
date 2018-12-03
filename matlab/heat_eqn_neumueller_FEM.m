%% solve heat eqn Ganda & Neumueller approach
% solve u_t - delta(u) = f
% parallelise implicit Euler, for delta FEM scheme
%%using weak formulation

% assemble system of L_h, M_h 

clear all;
close all;

tau = 0.04;
T = 2;
time_steps = ceil(T/tau)+1;
t = linspace(0,T, time_steps);

N = 20; 
h = 1/(N-1);
pts = N;

x = linspace(0,1, pts);

L_h = set_up_L_h_FEM(N);
M_h = set_up_M_h_FEM(N);

u0 = sin(2*pi*x)+1;

% set up f
f = zeros(pts, time_steps);
for j=1:time_steps
    for i=1:pts
    f(i,j) = 10*(sin(2*pi*x(i)) + cos(2*pi*t(j)));
    end
end

% add in boundary conditions
left_bdy = ones(time_steps, 1);
right_bdy = ones(time_steps, 1);

rhs = tau* M_h * f;

rhs(1:N, 1) = u0;
rhs(1,:) = left_bdy;
rhs(end,:) = right_bdy;

rhs_vec = reshape(rhs, [(pts*time_steps),1]);

% create list of blks
A_blks = {};
B_blks = {};

%% set L_h(1,1), L_h(end,end) to zero, that we get 1 when we add them
L_h(1,1) = 0; L_h(end,end) = 0;

%% assembly overall matrix
for ts=1:time_steps
    A_blks{ts} = M_h + tau*L_h;
    B_blks{ts} = -M_h;
    B_blks{ts}(1,1) = 0;
    B_blks{ts}(end,end) = 0;
end

A_blks{1} = eye(pts);

% assemble matrices and vectors
A_large = blkdiag(A_blks{:});
B_large = blkdiag(B_blks{1:end-1});
B_large = [zeros(pts, (time_steps-1)*pts); B_large];
B_large = horzcat(B_large, zeros(time_steps*pts, pts));

tot_mat = A_large + B_large;

%spy(tot_mat)
%u_large = tot_mat \ rhs_vec;

% not symmetric! need different method ...
% MULTIGRID PARAMETERS
u_init = [u0'; zeros((time_steps-1)*pts, 1)];
max_iter = 50;
smoother = 'GaussSeidel';

u_large = V_cycle(tot_mat, rhs_vec, u_init, 3, max_iter, smoother);

u_mat = zeros(pts, time_steps);

for ts=1:time_steps
    u_mat(:, ts) = u_large((ts-1)*pts+1:ts*pts);
end

figure
mesh(x, t, u_mat')
xlabel('space');
ylabel('time');
zlabel('u');
title('heat equation Neumueller FEM');




