%% solve heat eqn Ganda & Neumueller approach
% solve u_t - delta(u) = f
% parallelise implicit Euler, for delta FD scheme
%% not using weak formulation

% assemble system of L_h, M_h 

% include boundary conditions?!

clear all;
close all;

tau = 0.02;
T = 2;
time_steps = ceil(T/tau)+1;
t = linspace(0,T, time_steps);

N = 30; 
h = 1/N;
pts = N+1;

x = linspace(0,1, pts);
temp = [0; ones(pts-2, 1); 0];
%temp = ones(pts,1);

L_h = 1/tau*spdiags(temp, 0, pts, pts) + set_up_Lh_1D(N, h);
%M_h = set_up_M_h_1D(N);
%B = 1/tau * spdiags([-temp, temp], -1:0, pts, pts);
off_diag = -1/tau*spdiags(temp, 0, pts, pts);
off_diag(1,1) = 0;
off_diag(end,end) = 0;

u0 = sin(2*pi*x)+1;

% make f so that it can change over time

f = zeros(time_steps, pts);
for j=1:time_steps
    for i=1:pts
    f(j,i) = 10*(sin(2*pi*x(i)) + cos(2*pi*t(j)));
    end
end

% add in boundary conditions
left_bdy = ones(time_steps, 1);
right_bdy = ones(time_steps, 1);

f(:,1) = left_bdy;
f(:, end) = right_bdy;

f(1,1:pts) = u0;
f_vec = reshape(f', [(pts*time_steps),1]);

%only consider inner part
%L_h = L_h(2:end-1, 2:end-1);
%M_h = M_h(2:end-1, 2:end-1);

% create list of blks
A_blks = {};
B_blks = {};


%% assembly overall matrix
for ts=1:time_steps
    %A_blks{ts} = B_h + L_h;
    A_blks{ts} = L_h;
    %B_blks{ts} = -B_h;
    B_blks{ts} = off_diag;
end

A_blks{1} = eye(pts);

% assemble matrices and vectors
A_large = blkdiag(A_blks{:});
B_large = blkdiag(B_blks{1:end-1});
B_large = [zeros(pts, (time_steps-1)*pts); B_large];
B_large = horzcat(B_large, zeros(time_steps*pts, pts));

tot_mat = A_large + B_large;
% u_init = zeros((time_steps+1)*pts,1);
% max_iter = 50;
% smoother = 'GaussSeidel';

%spy(tot_mat)
u_large = tot_mat \ f_vec;
%[u, ~,~] = W_cycle(tot_mat, f_vec, u_init, 3, max_iter, smoother);

u_mat = zeros(time_steps, pts);
for ts=1:time_steps
    u_mat(ts, :) = u_large((ts-1)*pts+1:ts*pts);
end

figure
mesh(x, t, u_mat)
xlabel('space');
ylabel('time');
zlabel('u');




