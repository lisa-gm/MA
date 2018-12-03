%% set up Wave Eqn space time parallel
% Gabis approach

clear all;
close all;

c = 0.5;
% assume Nx = Nh, dx = dt
N = 10;

h = 1/(N-1);

x = linspace(0,1,N);
%x = x(2:end-1);

t = linspace(0,1,N);
%t = t(2:end-1);

total_pts = N^2;
s = c^2;

u0 = sin(2*pi*x);
%% set up K_h matrix
% walk first in time, then in space
temp = ones(total_pts,1);
I = eye(N);
Mid = spdiags([-s*temp, 2-2*s*temp, -s*temp], -1:1, N, N);
% fix first and last value through bdy conditions
Mid(1,:) = [1, zeros(1,N-1)];
Mid(end,:) = [zeros(1,N-1), 1];

Sides = spdiags([temp, temp], [-1, 1], N, N);
K_h = (kron(I, Mid) + kron(Sides, I));

K_h(1:N, :) = zeros(N, total_pts);
K_h(end-N+1:end, :) = zeros(N, total_pts);

K_h(1:N, 1:N) = I;
K_h(end-N+1:end, end-N+1:end) = I;

%% set up rhs

f = zeros(total_pts,1);
f(1:N) = u0;
f(end-N+1:end) = -u0;

%% SOLVE

u = K_h \ f;

norm(K_h*u - f)
%% plot solution

u_mat = reshape(u, [N, N]);
u_mat = u_mat';

figure
mesh(x, t, u_mat);
xlabel('space');
ylabel('time');
zlabel('u');

