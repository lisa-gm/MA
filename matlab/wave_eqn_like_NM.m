%% new try wave eqn parallelised 
% this time using approach like Neumueller
% except that I'm using a 2nd order backward difference
% to avoid problem that we need boundary conditions at the end

% what I want to solve: u_ij where i is space, j time
% 1/dt^2*(u_ij - 2u_ij-1 + u_ij-2) - c^2/dx^2 (u_i+1j-1 -2uij-1 + ui-1j-1)
% = f_ij

% for entry i,j we have 1/dt^2 - 2*c/dx^2, of diagonals -c/dx^2
% then two blocks in lower left
% one with 1/dt^2 other one with -2/dt^2

% include boundary conditions in everything

clear all;
% close all;

%% set up parameters
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

% bdy conditions, also make them change over time
bdy_l = zeros(Nt, 1);
bdy_r = zeros(Nt, 1);

%bdy_l = sin(2*pi*x);
%bdy_r = sin(2*pi*x);

% add in f, for now time independent
%f = 0.02*cos(pi*x);
%f = zeros(Nx,1);

u0 = sin(pi*x);
%u0 = sin(8*pi*x);
%u0(14:end) = 0;

%% set up the matrix
%K_h = set_up_K_h_2nd_or_backDiff(Nt, Nx, c2, ht, hx);
K_h = set_up_K_h_backDiff(Nt, Nx, c2, ht, hx);

% assemble rhs
rhs = zeros(Nx*Nt,1);

% put in initial conditions
% analytical solution is sin(pi*x)*cos(2*pi*t)

%% we need initial conditions for t0 and t1
rhs(1:Nx) = u0;
rhs(Nx+1:2*Nx) = u0;
%rhs(Nx+1:2*Nx) = sin(pi*x)*cos(2*pi*t(2));
%rhs(2*Nx+1:3*Nx) = sin(pi*x)*cos(2*pi*t(3));

% put in boundary conditions
rhs(1:Nx:end) = bdy_l;
rhs(Nx:Nx:end) = bdy_r;

%% now solve 

% would need a solver that can handle non-symmetric matrices
u = K_h \ rhs;

%% get result into matrix format

u_mat = reshape(u, [Nx, Nt]);
u_mat = u_mat';

res = rhs - K_h*u;
err = zeros(1, Nt);
for ts = 1: Nt
    err(ts) = norm(res((ts-1)*Nx+1:ts*Nx));
end

sum(err)

%% plot result in one
figure
mesh(x, t, u_mat);
xlabel('space');
ylabel('time');
zlabel('u');
title('wave eqn like Neumueller');

%% plot result as time lapse
% 
% figure
% plot(x, u_mat(1,:));
% axis([0 1 -1 2]);
% 
% for ts=1:Nt
%     plot(x, u_mat(ts, :));
%     axis([0 1 -1 2]);
%     pause(0.1);
% end



