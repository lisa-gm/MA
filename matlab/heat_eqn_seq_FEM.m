%%%%% incorporate FE solver and implicit Euler %%%%%%%%%
%%% to solve du/dt - delta(u) = f

clear all;
close all;

Nx = 100; 
K = 1;

h = 1/(Nx-1);
x = linspace(0,1, Nx);

T = 2;
dt = 10^(-1);
Nt = T/dt;

t = linspace(0, T, Nt);

% multigrid parameters
max_iter = 5000;
levels = 3;

% u(0,x):
%u0 = ones(Nx, 1);
u0 = sin(pi*x)+1;
%u0 = transpose(1/(4*pi^2) * sin(2*pi*x)+1);
%u0 = transpose(2* sin(pi*x/2) - sin(pi*x) + 4*sin(2*pi*x));

% solve for specific time t in space with FE
% set up discretised Laplacian

L_h = set_up_L_h_FEM(Nx);
M_h = set_up_M_h_FEM(Nx);

%%%% if we have additional term f %%%%%
%f = transpose(-10*(x-0.5).^2+0.25);
%f = 10*ones(length(x),1);
f = zeros(Nx, Nt);
for j=1:Nt
    f(:,j) = 0; %sin(pi*x)cos(pi*t(j));
    %f(i,j) = 10*(sin(2*pi*x(i)) + cos(2*pi*t(j)));
end

% add in boundary conditions
left_bdy = ones(Nt, 1);
right_bdy = ones(Nt, 1);

f(1,:) = left_bdy;
f(end, :) = right_bdy;

f(1:Nx,1) = u0;

f_vec = reshape(f, [(Nx*Nt),1]);

% Neumann boundary conditions, du/dn = ... on dOmega
% need initial value at beginning, otherwise off by constant
%NM = [0, 0];
%f(1)= f(1) + NM(1);
%f(end) = f(end) + NM(2);

% if we use multigrid instead of computing inverse directly
% want to solve (M+dt*L_h)*u_t+1 = M * u_t + delta_t*M_h*f

y = zeros(Nx, Nt);
y(:,1) = u0;

% for IMP Euler
A = M_h+dt*L_h;

% for Crank Nicolson
%k = 0.0001;
%A = M_h + k/2*L_h;
 
% set bc
A(1,1) = 1;
A(end,end)= 1;

for st=2:Nt
    % for IMP EULER    
    rhs =  M_h*y(:,st-1) + dt*M_h*f(:, st); 
    
    % for Crank Nicolson
    %k = 0.001;
    %rhs = (M_h - k/2 * A)*y(:,st-1);
    
    rhs(1) = left_bdy(st);
    rhs(end) = right_bdy(st);
    
    A_Gauss = A;
    
    rhs(2) = rhs(2) - A(2,1)*rhs(1);
    rhs(end-1) = rhs(end-1) - A_Gauss(end-1,end)*rhs(end);
    A_Gauss(2,1) = 0;
    A_Gauss(end-1,end) = 0;   
    
    y(:, st) = A_Gauss \ rhs;
    
    %[y(:,st)] = V_cycle(A_Gauss, rhs, y(:,st-1), levels, max_iter, 'GaussSeidel');
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
title('heat equation seq FEM');
%plot(x, y(:,1), x, y(:,floor(no_time_st/4)), x, y(:,floor(no_time_st/2)), x, y(:,3*floor(no_time_st/4)), x, y(:,end));
%legend('u0', 'u1', 'u2', 'u3', 'u end');

% analytical solution
% u_exact = zeros(Nx, Nt);
% 
% for ts=1:Nt
%     u_exact(:, ts) = sin(pi*x)*exp(-pi^2*t(ts))+1;
% end
% 
% norm(u_exact - y)
% 
% figure
% mesh(x, t, u_exact');
% xlabel('space');
% ylabel('time');
% zlabel('u');
% title('heat equation analytical');

