% *********************************************************** %
% ************************ RUN LSFEM ************************ %
% *********************************************************** %

clear all;
close all;

% CHOOSE GENERAL PARAMTERS: 

max_iter = 30;                          % number of newton iterations
omega = 1;                              % for a scaled newton step
eps = 10^-12;                           % tolerance

alpha = 1;                             % scaling in front of f(u) ie mu*f(u) -> derivatives
diff_const = 10^(-3);                   % diffusion constant, define through sigma = diff_const*nabla(u) 
                                                                               
c1 = 2;                                 % constants in front of J = c1*|| ....-f(u) ||^2 + c2*|| ...||^2
c2 = 2;                                 % need to multiply with them in computation of each integral 
                                        
% *********** DOMAIN SET UP *************** 
S = 1;
Nx_elem = 5;

T = 1;
Nt_elem = 5;

% SET UP DOMAIN, to compute the rest ... 
int_space = [0,S];
hx = (int_space(2) - int_space(1))/Nx_elem;
Nx_pts = Nx_elem+1;
x_vec = linspace(int_space(1), int_space(2), Nx_pts);

int_time = [0, T];
ht = (int_time(2) - int_time(1))/Nt_elem;
Nt_pts = Nt_elem+1;
t_vec = linspace(int_time(1), int_time(2), Nt_pts);

tot_elem = Nx_elem*Nt_elem;
tot_pts = Nx_pts*Nt_pts;

hxht = ht*hx;

% ********* COMPUTING EXACT SOLUTION (if possible) **********

u_exact = zeros(Nx_pts, Nt_pts);
sigma_exact = zeros(Nx_pts, Nt_pts);

for ts=1:Nt_pts
    %u_exact(:,ts) = cos(pi*t(ts))*sin(pi*x); % x.^2*t(ts); 
    %u_exact(:, ts) = 5;
    %u_exact(:, ts) = sin(2*pi*t_vec(ts)); %*cos(pi*t_vec(ts));
    %u_exact(:, ts) = t_vec(ts);
    %u_exact(:,ts) = x_vec;
    %u_exact(:,ts) = t_vec(ts)*(x_vec).^2;
    
    %sigma_exact(:, ts) = diff_const*2*pi*cos(2*pi*x_vec);
    %sigma_exact(:,ts) = diff_const*2*t_vec(ts)*x_vec;
    %sigma_exact(:, ts) = diff_const*pi*cos(pi*x_vec)*cos(pi*t_vec(ts));
    %sigma_exact(:, ts) = diff_const*1;
end

% *************** BOUNDARY CONDITIONS ******************
%bdy_cond = 'Dirichlet';
bdy_cond = 'Neumann';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(strcmp(bdy_cond, 'Dirichlet'))
    %bdy_left = sin(2*pi*t_vec);
    bdy_left = 0*ones(Nt_pts, 1);
    %bdy_left = t_vec;
    
    %bdy_right = sin(2*pi*t_vec);
    %bdy_right = 1*ones(Nt_pts, 1);
    bdy_right = t_vec;
end

% assume neumann bdy cond of the form \sigma*n_l = g_l, n_l = -1
% \sigma*n_r = g_r, n_r = 1
% enforce this in the matrix
if(strcmp(bdy_cond, 'Neumann'))
    bdy_left = zeros(Nt_pts,1);
    bdy_right = zeros(Nt_pts,1);
    
    for ts=1:Nt_pts
        bdy_left(ts) = 0;
        bdy_right(ts) = 0;
    end  
end

% ************** initial conditions *****************
%u0 = x_vec.^2;
%u0 = sin(2*pi*x_vec);
%u0 = x_vec;
%u0 = 0;
u0 = 1 - 2.*x_vec;
u0(floor(Nx_pts/2)+1:end) = 0;
%u0 = u_exact(:,1);

u = ones(tot_pts, 1);

u(1:Nx_pts) = u0;
u(1:Nx_pts:end) = bdy_left;
u(Nx_pts:Nx_pts:end) = bdy_right;

sigma = ones(tot_pts,1);

sol = [sigma; u];

% ********** define f ***************
% polynomial of the form: sum_i a_i u^i with a_i = a(x,t), i<=3
% f = a0 + a1*u + a2*u^2 + a3*u^3
% a0 = @(x,t) 0; a1 = @(x,t) 0; a2 = @(x,t) 0; a3 = 0;

%a0 = @(x,t) - pi*sin(pi*x)*sin(pi*t) + pi^2*sin(pi*x)*cos(pi*t); %- (sin(pi*x)*cos(pi*t))^3;
%a0 = @(x,t) x*x - diff_const*2*t - (x^2*t)^2;
%a0 = @(x,t) - x^2;
%a1 = @(x,t) 0; 
%a2 = @(x,t) 1; 
%a3 = 0; 

 a0 = @(x,t) 0;
 a1 = @(x,t) 0; %-0.1; 
 a2 = @(x,t) 0; %1.1; 
 a3 = 0; %1;

% ***************** call functional J ********************** %
%sigma = reshape(sigma_exact, [tot_pts, 1]);
%u = reshape(u_exact, [tot_pts, 1]);

%J_eval = J_eval_lsfem(S, Nx_elem, T, Nt_elem, sigma, u, a0, a1, a2, a3, alpha, bdy_cond, c1, c2, diff_const);

% ********************* CALL NONLINEAR SOLVER *******************

[sigma_mat, u_mat] = lsfem_nonlinear_fct(max_iter, omega, eps, alpha, diff_const, S, Nx_elem, T, Nt_elem, ...
    u_exact, sigma_exact, bdy_cond, bdy_left, bdy_right, u0, a0, a1, a2, a3, c1, c2);