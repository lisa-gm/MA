%%%%% incorporate FE solver and implicit Euler %%%%%%%%%
%%% to solve delta_x = du/dt 

clear all;
close all;

Nx = 50; 
K = 1;

hx = 1/(Nx-1);
x = linspace(0,1, Nx);

T = 5;
dt_coarse = 0.5;
Nt = T/dt_coarse;

t = linspace(0, T, Nt);

% now refine in each step, no of fine steps in each interval
no_time_st_fine = 10;
dt_fine = dt_coarse/no_time_st_fine;

% multigrid parameters
max_iter = 30;
levels = 3;

% u(0,x):
%u0 = ones(Nx, 1);
u0 = transpose(1/(4*pi^2) * sin(2*pi*x)+1);
%u0 = transpose(2* sin(pi*x/2) - sin(pi*x) + 4*sin(2*pi*x));
% check if they coincide with boundary conditions over time

% solve for specific time t in space with FE
% set up discretised Laplacian
L_h = set_up_L_h_FEM(Nx);
M_h = set_up_M_h_FEM(Nx);

%%%% if we have additional term f %%%%%
%f = transpose(-(x-0.5).^2+0.25);
%f = 10*ones(length(x),1);
f = zeros(Nx, Nt);
%f(ceil(N/2)) = 4;

% for Dirichlet boundary conditions
for j=2:Nt
    for i=1:Nx
    f(i,j) = 10*(sin(2*pi*x(i)) + cos(2*pi*t(j)));
    end
end

% Neumann boundary conditions, du/dn = ... on dOmega
% need initial value at beginning, otherwise off by constant
%NM = [1, 1];
%f(1)= f(1) + NM(1);
%f(end) = f(end) + NM(2);

% if we use multigrid instead of computing inverse directly
% want to solve (M+delta_t*A)*u_t+1 = M * u_t + delta_t*M_h*f

y = zeros(Nx, Nt);
y(:,1) = u0;

% for IMP Euler
A = M_h+dt_coarse*L_h;

% for Crank Nicolson
%k = 0.0001;
%A = M_h + k/2*L_h;
 
% set bc
A(1,1) = 1;
A(end,end)= 1;

%%%%% OUTER LOOP FOR OVERALL ITERATION %%%%%%


% if i want to parallelise in time by first computing 
% a few big time steps and then use those for 
% simultaneous refinement

for st=2:Nt
    % for IMP EULER
    rhs =  M_h*y(:,st-1) + (dt_coarse*M_h*f(:,st)); 
    
    % for Crank Nicolson
    %rhs = (M_h - k/2 * A)*y(:,st-1);
    
    rhs(1) = u0(1);
    rhs(end) = u0(end);
    
    A_Gauss = A;
    
    rhs(2) = rhs(2) - A(2,1)*rhs(1);
    rhs(end-1) = rhs(end-1) - A_Gauss(end-1,end)*rhs(end);
    A_Gauss(2,1) = 0;
    A_Gauss(end-1,end) = 0;   
    
    y(:, st) = A_Gauss \ rhs;
    
    %[y(:,st), ~, ~] = W_cycle(A_Gauss, rhs, y(:,st-1), levels, max_iter, 'GaussSeidel');
end

% now we have the coarse approximation, use for refinement
% use y(:, st) as initial values

% rewrite everything into matrices to perform actions in parallel

rhs_fine = zeros(Nx, no_time_st_fine);
y_fine = zeros(Nx, no_time_st_fine*(Nt-1) + 1);

for st_fine=2:no_time_st_fine
        % for IMP EULER
    rhs_fine =  M_h*y + dt_fine*M_h*f; 
    
    % for Crank Nicolson
    %rhs = (M_h - k/2 * A)*y(:,st-1);
    
    rhs_fine(1,:) = u0(1);
    rhs_fine(end,:) = u0(end);
    
    A_Gauss = A;
    
    rhs_fine(2,:) = rhs_fine(2,:) - A_Gauss(2,1)*rhs_fine(1,:);
    rhs_fine(end-1,:) = rhs_fine(end-1,:) - A_Gauss(end-1,end)*rhs_fine(end,:);
    A_Gauss(2,1) = 0;
    A_Gauss(end-1,end) = 0;   
    
    %y(:, st) = A_Gauss \ rhs;
    
    % multigrid can't handle matrix input, just put loop now
    % however this can be done in parallel, independent of each other
    for i=1:Nt-1
    [y_fine(:,no_time_st_fine*(i-1)+st_fine), ~, ~] = V_cycle(A_Gauss, rhs_fine(:,i), y_fine(:,no_time_st_fine*(i-1)+st_fine-1), levels, max_iter, 'GaussSeidel');
    end
end

% we now have old and new entries for the coarse level steps
% just average them out

for i=1:(Nt-1)
y_fine(:,i*no_time_st_fine+1) = (y_fine(:,i*no_time_st_fine+1) + y(:,i+1))./2;
end

t_fine = linspace(0,T, no_time_st_fine*T);

size(y_fine)
size(t_fine)

%%% try with Crank Nicolson 
%k = 0.001;
%y = CN_method(M_h, L_h, u0, no_time_st, k);

u_t = y(:,end);

figure;
%plot(x, y(:,1), x, y(:,floor(no_time_st/4)), x, y(:,floor(no_time_st/2)), x, y(:,3*floor(no_time_st/4)), x, y(:,end));
%legend('u0', 'u1', 'u2', 'u3', 'u end');
%plot(x, y_fine(:,end), x, y(:,end), x, y(:,1))
% xlim([0 1]);
% ylim([0 1]);
% legend('y fine', 'y coarse end', 'y coarse initial');
mesh(x,t_fine, y_fine')
xlabel('space');
ylabel('time');
zlabel('u');
title('heat equation parareal');





