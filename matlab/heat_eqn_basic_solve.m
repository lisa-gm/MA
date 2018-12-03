%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% SETTING UP LEAST SQUARES SPACE TIME HEAT EQ %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%% u \in P1P1 AND \sigma in P0P1 %%%%%%%%%%%%%%%%%%%%%
%clear all;
%close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% GENERAL PARAMETERS 

Nx = 20;
S = 1;
hx = S/(Nx-1);
x = linspace(0, S, Nx);

Nt = 80;
T = 5;
ht = T/(Nt-1);
t = linspace(0, T, Nt);

%u0 = transpose(-2.5*x.^2 + 2.5*x+1);
u0 = transpose(1/(pi^2)*sin(pi*x)+1);
%u0 = zeros(Nx, 1);

% f
%f = zeros(Nx, Nt);
f = ones(Nx, Nt);

for ts=1:Nt
    %f(:,ts) = 5;
    f(:,ts) = sin(pi*x)*cos(pi*t(ts)); %+ cos(2*pi*t(ts));
   
end

f_vec = reshape(f, [(Nx*Nt),1]);

%%%%%%%%% BOUNDARY CONDITIONS

bdy_left = ones(Nt, 1);
bdy_right = ones(Nt, 1);
%bdy_left = zeros(Nt, 1);
%bdy_right = zeros(Nt, 1);
%bdy_left = linspace(0,1, Nt);
%bdy_right = linspace(0,1,Nt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% GET STIFFNESS, GRADIENT, MASS Matrices in time + space

[G_h, K_h, M_h] = set_up_FEM_1D_mat(Nx);
[G_t_u, K_t_u, M_t_u] = set_up_FEM_1D_mat(Nt);
[G_p0p1, M_p0p1, M_p0p0] = set_up_P0P1_FEM_mat(Nt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% SOLVE by using FEM in space, implicit Euler in time
%%%%%%% WRITE AS BIG SYSTEM

% A_block_seq = M_h + ht*K_h;
% 
% temp = ones(Nt,1);
% A_seq = kron(eye(Nt), A_block_seq) + kron(spdiags(-1*temp, -1, Nt, Nt), M_h);
% 
% A_seq(1:Nx, :) = 0; A_seq(1:Nx:end, :) = 0; A_seq(Nx:Nx:end) = 0;
% 
% A_seq(1:Nx, 1:Nx) = eye(Nx);
% A_seq(1:Nx:end, 1:Nx:end) = eye(Nt); 
% A_seq(Nx:Nx:end, Nx:Nx:end) = eye(Nt);
% 
% f_int = ht*M_h*f;
% f_int(:,1) = u0;
% f_int(1,:) = bdy_left;
% f_int(end,:) = bdy_right; 
% 
% f_int_vec = reshape(f_int, [Nx*Nt,1]);
% 
% u_seq = A_seq \ f_int_vec;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%  PRIMAL SPACE TIME FORMULATION FEM %%%%%%%%%%%%%%%%%%%%%%%%%%%%

area = T*S;
%scale = 1/(ht*hx);
scale = 1;
A_primal = full(scale*kron(G_t_u, M_h) + area*scale*kron(M_t_u, K_h)); %


A_primal(1:Nx,:) = 0;
A_primal(1:Nx:end, :) = 0;
A_primal(Nx:Nx:end, :) = 0;

A_primal(1:Nx, 1:Nx) = eye(Nx);
A_primal(1:Nx:end, 1:Nx:end) = eye(Nt);
A_primal(Nx:Nx:end, Nx:Nx:end) = eye(Nt);

f_primal = area*scale*kron(M_t_u, M_h)*f_vec;

f_primal(1:Nx) = u0;
f_primal(1:Nx:end) = bdy_left;
f_primal(Nx:Nx:end) = bdy_right;

u_primal = A_primal\f_primal;


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%% PLOT

u_mat_seq = zeros(Nx, Nt);
u_mat_prim = zeros(Nx, Nt);


for ts=1:Nt
    %u_mat_seq(:, ts) = u_seq((ts-1)*Nx+1:ts*Nx);
    u_mat_prim(:, ts) = u_primal((ts-1)*Nx+1:ts*Nx);
end

figure
mesh(x, t(1:end), transpose(u_mat_prim(:,1:end)));
xlabel('space');
ylabel('time');
zlabel('u');
title('heat equation LSFEM');
% 
% figure
% 
% subplot(1,2,1);
% mesh(x, t(1:end), transpose(u_mat_seq(:,1:end)));
% xlabel('space');
% ylabel('time');
% zlabel('u');
% title('heat equation sequential formulation');

% subplot(1,2,2);
% mesh(x, t(1:end), transpose(u_mat_prim(:,1:end)));
% xlabel('space');
% ylabel('time');
% zlabel('u');
% title('heat equation primal space-time formulation');

% f_exp = 2.*t';
% %lambda = -0.5; 
% f_exp_int = M_t_u*f_exp;
% 
% A_exp = G_t_u;
% A_exp(1,:) = 0; 
% A_exp(1,1) = 1;
% 
% f_exp(1,1) = 0;
% 
% u_exp = A_exp \ f_exp;
% 
% figure
% plot(t, u_exp, t, t.^2);
% legend('approx sol', 'real sol');
% title('variational formulation deriv(u) = lambda*u');
% xlabel('time');


