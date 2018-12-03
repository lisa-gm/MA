%%%%% incorporate FE solver and implicit Euler %%%%%%%%%
%%% to solve du/dt - delta(x) = f 

clear all;
close all;

Nx = 50; 
K = 1;

hx = 1/(Nx-1);
x = linspace(0,1, Nx);

T = 5;
dt = 10^(-2);
Nt = T/dt;

t = linspace(0, T, Nt);

% multigrid parameters
max_iter = 500;
levels = 3;

% u(0,x):
%u0 = ones(Nx, 1);
u0 = 1/(2*pi)^2*sin(2*pi*x)+1;
%u0 = transpose(1/(4*pi^2) * sin(2*pi*x)+1);
%u0 = transpose(2* sin(pi*x/2) - sin(pi*x) + 4*sin(2*pi*x));

L_h = set_up_L_h_FD(Nx);

%%%% if we have additional term f %%%%%
%f = transpose(-10*(x-0.5).^2+0.25);
%f = 10*ones(length(x),1);
f = zeros(Nx, Nt);
for j=2:Nt
    for i=1:Nx
    f(i,j) = (sin(2*pi*x(i)) + cos(2*pi*t(j)));
    end
end

% add in boundary conditions
left_bdy = ones(Nt, 1);
right_bdy = ones(Nt, 1);

f(1, :) = left_bdy;
f(end, :) = right_bdy;

f_vec = reshape(f, [Nx*Nt,1]);

% Neumann boundary conditions, du/dn = ... on dOmega
% need initial value at beginning, otherwise off by constant
%NM = [0, 0];
%f(1)= f(1) + NM(1);
%f(end) = f(end) + NM(2);

y = zeros(Nx, Nt);
y(:, 1) = u0;

M = eye(Nx, Nx) + dt*L_h;

for st=2:Nt
    % IMP EULER
    
    % for Crank Nicolson
    %rhs = (M_h - k/2 * A)*y(:,st-1);
    
    rhs = y(:,st-1) + dt*f(:, st);
    M_Gauss = M;
   
    % to make matrix symmetric
    rhs(2) = rhs(2) - M(2,1)*rhs(1);
    rhs(end-1) = rhs(end-1) - M_Gauss(end-1,end)*rhs(end);
    M_Gauss(2,1) = 0;
    M_Gauss(end-1,end) = 0;   
    
    y(:,st) =  M \ rhs; 
    
    %[y(:,st), ~, ~] = W_cycle(M_Gauss, rhs, y(:,st-1), levels, max_iter, 'GaussSeidel');
end

%%% try with Crank Nicolson 
%k = 0.001;
%y = CN_method(M_h, L_h, u0, no_time_st, k);

u_t = y(:,end);

figure
mesh(x, t, y');
xlabel('space');
ylabel('time');
zlabel('u');
title('heat equation seq FD');
