%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% 1D MONODOMAIN EQUATION %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clear all;


% taken from Adaptive numerical sol. for PDEs, Deuflhard&Weiser
% du/dt = c*d2u/dx2 + alpha*f(u)
% f(u) = u*(1-u)*(u-0.1) = -u^3 + 1.1*u^2 - 0.1*u
% f  : IR -> IR
% f' : IR -> IR

% solve with space-time-matrix + weak formulation
% rewrite as A*u = alpha*f(u) 
% TO DO: integrate c into matrix ...

% solve with Newton, 
% need initial iterate u_guess on entire space-time domain
% J(u) = alpha*f(u) - A*u = 0
% u is vector -> pointwise evaluations of f, gives vector
% H(u) = grad(g(u)) = diag(alpha*f'(u)) - A

%% u_k+1 = u_k - (grad(g(u)k)))^-1)*g(u_k)


%%%%%%%%%%%%%% PARAMETER INITIALISATION %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% EQUATION PARAMETERS

% c = 10^-3;
% a = 10;

c = 10^-3;
%c = 1;
a = 10; 

% NEWTON SOLVER PARAMETERS
iter = 1;
max_iter = 20;
eps = 10^-9;
omega = 1;

Nx = 21;
hx = 1/(Nx-1);
x = 0:hx:1;

T = 2;
Nt = 21;
ht = T/(Nt-1);
t = 0:ht:T;

% initial conditions at t=0
u0 = zeros(Nx,1);
for i=1:Nx
    u0(i) = max(0, 1-2*x(i));
end
% u0 = 5*ones(Nx,1);

% how to evaluate f efficiently?
f = @(u) u.*(1-u).*(u-0.1);
deriv_f = @(u) -3*u.*u + 2.2*u - 0.1;

% change f a bit so something happens
% f = @(u) u.*(2-u).*(u-0.1);             % = (2u^2 - u^3 - 0.2u  )(u - 0.1)
% deriv_f = @(u) -3*u.*u + 4.2*u - 0.2;   

%% LINEAR TESTCASE
% f = @(u) [u0; zeros(length(u)-length(u0), 1)];
% deriv_f = @(u) zeros(length(u),1);

% homogenous Neumann == No boundary conditions !
% bdy_left = zeros(Nt, 1);
% bdy_right = zeros(Nt, 1);
% 
% u0 = [bdy_left(1); u0(2:end-1); bdy_right(1)];

%%%%%%%%%%%%%%%%%%%%% SET UP MATRICES %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[comp_mat, K_h, K_t, M_h, M_t, N_t] = set_up_M_st(ht, Nx, Nt, c, 'Neumann');

%%%%%%%%%%%%%%%%%%%% NEWTON SOLVER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialisation u_guess
N_tot = size(comp_mat, 1);
u = 0.1*ones(N_tot, 1);
u(1:Nx) = u0;

% bdy_left_DG = repelem(bdy_left, 2);
% bdy_left_DG = bdy_left_DG(2:end-1);
% bdy_right_DG = repelem(bdy_right,2);
% bdy_right_DG = bdy_right_DG(2:end-1);
% 
% u(1:Nx:end) = bdy_left_DG;
% u(Nx:Nx:end) = bdy_right_DG;

%% while loop to solve
% just initialise some M_f so that we go into while loop
M_f = 10*ones(length(u),1);

%%%%%%%%%%%%%%%%%%%%%%%% LOOP %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% another stopping criterion: norm(u_new - u_old) < eps
u_old = -1*ones(length(u),1);
u_direct_old = u_old;

%%% MG parameters:
levels = 3;
smoother = 'GMRES';
max_iter_mg = 4;

e = 0.1*ones(size(u));

while(iter <= max_iter && norm(u - u_old) > eps)
    tic;
    f_mat = reshape(f(u), [2*Nx, Nt-1]);
    M_f_mat = kron(M_t, M_h)*f_mat;
    M_f = reshape(M_f_mat, [2*Nx*(Nt-1),1]);
    
    M_f(1:Nx) = M_h*u0;
    
    J = comp_mat*u - a*M_f;

    % just define first entries to be zero, seems to work better
    J(1:Nx) = zeros(Nx, 1);
    
    deriv_f_mat = reshape(deriv_f(u),[2*Nx, Nt-1]);
    M_deriv_f_mat = kron(M_t, M_h)*deriv_f_mat;
    M_deriv_f = reshape(M_deriv_f_mat, [2*Nx*(Nt-1),1]);
    
    % also IC business in derivative, at least for NL case?
    H = sparse(comp_mat - diag(a*M_deriv_f));
    
    % check if H diagonally dominant
%     arg_max = zeros(1, size(H,1));
%     for row=1:size(H,1)
%         [~, arg_max(row)] = max(abs(H(row,:)));
%     end
%     
%     
%     temp = arg_max - 1:size(H,1);
%     temp(temp == 0) = 1;
%     temp(temp ~= 0) = 0;
    
    
    u_old = u;
    
    % direct solve
    % u_direct = u_old - omega*(H \ J);
    % u_direct_old = u_direct;

    %u = u_direct;
    
    e = H \ (-J);
    
    % only smoothing
    % smoothing(H, -J, e, 5, smoother)
    
    % geometric multigrid
    %e = V_cycle(H, -J, e, levels, max_iter_mg, smoother);

    
    % AMG
    %e = two_level_AMG(H, -J, e, smoother, eps, max_iter_mg);

    u = u + e;

    iter = iter+1;
    tt = toc;
    fprintf('one complete Newton step incl AMG, time: %d sec, iter: %d \n', tt, iter-1);
end


%% for linear case, can also use direct solve to compare

% u0_rhs = kron([1,0],u0');  
% M_f(1:2*Nx) = kron([1 0 ; 0 0],M_h) * u0_rhs';
% 
% u = comp_mat \ M_f;

%% Analysis of the results
fprintf('number of iterations: %d\n', iter-1);

% find better way to compute error
fprintf('norm residual: %d\n', norm(J, inf));

%% rearrange
% get 2 values for each time step 
% better to take the first one ?!
u_mat = zeros(Nx, Nt);
% u_mat_direct = zeros(Nx, Nt);

u_mat(:,1) = u0;
%u_mat_direct(:,1) = u0;
% note that we are only taking one set of solutions for 
% the inner nodes, always x_i- 
% WHY better than x_i+ ?? 

% indices x_i+
% u_mat(:, ts) = u((2*ts-2)*Nx+1:2*(ts-1)*Nx+Nx)?

for ts=2:Nt
    %u_mat_direct(:,ts) = u_direct((2*ts-3)*Nx+1:2*(ts-1)*Nx);
    u_mat(:, ts) = u((2*ts-3)*Nx+1:2*(ts-1)*Nx);
end
%save('u_approx_sol.mat', 'u');
%save('H.mat', 'H', 'u_mat','x','t', 'Nx', 'Nt', 'N_tot');

%%%%%%%%%%%%%%%% VISUALIZING THE RESULTS %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
suptitle('approximated solution');

%subplot(1,2, 1);
mesh(x, t, u_mat'); 
xlabel('space');
ylabel('time');
% subplot(2,2, 2);
% mesh(x, t, u_mat_direct'); 
% xlabel('space');
% ylabel('time');
