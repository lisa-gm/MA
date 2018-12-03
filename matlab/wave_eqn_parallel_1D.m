%% Wave Equation %% 
% consider d2u/dt2 = c^2*d2u/dx2
% rewrite in "laplacian form" and take finite difference scheme
% to build matrix 

close all;
clear all;

%% set up parameters

% for now same number of steps in space and time

c2 = 4;
c = sqrt(c2);

int_time = 10;
int_space = 1;

Nx = 50;
Nt = 500;

hx = int_space/(Nx-1);
ht = int_time/(Nt-1);

x = linspace(0,int_space,Nx);
t = linspace(0,int_time,Nt);

% bdy conditions
bdy_l = zeros(Nt, 1);
bdy_r = zeros(Nt, 1);

% add in f, for now time independent
%f = 0.02*cos(pi*x);
%f = zeros(Nx,1);

u0 = sin(pi*x);
u_end = sin(pi*x);
%u0(10:end) = 0;

K_h = set_up_K_h(Nt, Nx, c2, ht, hx);
%K_h = set_up_K_h_backDiff(N, p, delta_x, delta_t);
spy(K_h)

% change first and last entry of each block to take bdy conditons
% as well as intial conditions

%f(1:pts:end) = u_bdy_l;
%f(pts:pts:end) = u_bdy_r;

rhs = zeros(Nx*Nt,1);
rhs(1:Nx,1) = u0;
rhs(end-Nx+1:end,1) = u_end;

%% solve
% use direct solve for now
u = K_h \ rhs;
%[u, ~, ~] = W_cycle(K_h, f, u_init, levels, max_iter_MG, 'GaussSeidel');

u_mat = reshape(u, [Nx, Nt]);
u_mat = u_mat';

max_val = max(max(u_mat));
min_val = min(min(u_mat));

%% plot
u_final_t = u_mat(end,:);

figure
%plot(x,u(1:Nx+1))
mesh(x, t, u_mat);
xlabel('space');
ylabel('time');
zlabel('u');

norm(K_h*u - rhs)
