%%%%% incorporate FE solver and implicit Euler %%%%%%%%%
%%% to solve delta_x = du/dt 

% what do we need as input?
% initial conditions on u(0,x)
% and boundary values for u for all t

% what is our domain and how do we discretize it?
% for now 1D, [0,1]

clear all;
close all;

N = 100; 
K = 1;

h = 1/N;
x = linspace(0,1, N+1);
% take only inner grid points
%x = x(2:end-1);

delta_t = 10^(-2);
no_time_st = 100;

% multigrid parameters
max_iter = 5000;
levels = 3;

% u(0,x):
u0 = ones(N+1, 1);
%u0 = transpose(1/(4*pi^2) * sin(2*pi*x)+1);
%u0 = transpose(2* sin(pi*x/2) - sin(pi*x) + 4*sin(2*pi*x));
% check if they coincide with boundary conditions over time

% solve for specific time t in space with FE
% set up discretised Laplacian
L_h = set_up_L_h_1D(N, 1);
M_h = set_up_M_h_1D(N);

%%%% if we have additional term f %%%%%
%f = transpose(-10*(x-0.5).^2+0.25);
%f = 10*ones(length(x),1);
f = zeros(N+1, 1);
%f(ceil(N/2)) = 4;

% for Dirichlet boundary conditions
f(1) = u0(1);
f(end) = u0(end);

% Neumann boundary conditions, du/dn = ... on dOmega
% need initial value at beginning, otherwise off by constant
NM = [0, 0];
f(1)= f(1) + NM(1);
f(end) = f(end) + NM(2);

% if we use multigrid instead of computing inverse directly
% want to solve (M+delta_t*A)*u_t+1 = M * u_t + delta_t*M_h*f

y = zeros(size(u0,1), no_time_st+1);
y(:,1) = u0;

% for IMP Euler
A = M_h+delta_t*L_h;

% for Crank Nicolson
%k = 0.0001;
%A = M_h + k/2*L_h;
 
% set bc
A(1,1) = 1;
A(end,end)= 1;

for st=2:no_time_st+1
    % for IMP EULER
    rhs =  M_h*y(:,st-1) + (delta_t*M_h*f); 
    
    % for Crank Nicolson
    %rhs = (M_h - k/2 * A)*y(:,st-1);
    
    rhs(1) = u0(1);
    rhs(end) = u0(end);
    
    A_Gauss = A;
    
    rhs(2) = rhs(2) - A(2,1)*rhs(1);
    rhs(end-1) = rhs(end-1) - A_Gauss(end-1,end)*rhs(end);
    A_Gauss(2,1) = 0;
    A_Gauss(end-1,end) = 0;   
    
    %y(:, st) = A_Gauss \ rhs;
    
    [y(:,st), ~, ~] = W_cycle(A_Gauss, rhs, y(:,st-1), levels, max_iter, 'GaussSeidel');
end

%%% try with Crank Nicolson 
%k = 0.001;
%y = CN_method(M_h, L_h, u0, no_time_st, k);

u_t = y(:,end);

figure;
plot(x, y(:,1), x, y(:,floor(no_time_st/4)), x, y(:,floor(no_time_st/2)), x, y(:,3*floor(no_time_st/4)), x, y(:,end));
legend('u0', 'u1', 'u2', 'u3', 'u end');






