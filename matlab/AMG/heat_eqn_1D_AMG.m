%% solve heat eqn Ganda & Neumueller approach %%%%
%%%%%%%% WITH ALGEBRAIC MULTIGRID %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% solve u_t - delta(u) = f
% parallelise implicit Euler, for delta FEM scheme
%%using weak formulation
% NOT USING DG IN TIME

% assemble system of L_h, M_h 

clear all;
close all;

c = 1;

%%%%%%%%%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% TIME
T = 3;
Nt = 80;
time_steps = Nt;
ht = T/(Nt-1);
t = linspace(0,T, time_steps);

%%% SPACE 
Nx = 10; 
h = 1/(Nx-1);
pts = Nx;
x = linspace(0,1, pts);

bdy_type = 'Neumann';

L_h = set_up_L_h_FEM(Nx, bdy_type);
M_h = set_up_M_h_FEM(Nx, bdy_type);

u0 = sin(pi*x)+1;
%u0 = 5*ones(1, Nx);

% set up f
f = zeros(pts, time_steps);
for j=1:time_steps
    for i=1:pts
    f(i,j) = 10*(sin(2*pi*x(i)) + cos(2*pi*t(j)));
    end
end


% add in boundary conditions
left_bdy = ones(time_steps, 1);
right_bdy = ones(time_steps, 1);


% create list of blks
A_blks = {};
B_blks = {};

%% set L_h(1,1), L_h(end,end) to zero, that we get 1 when we add them
if(strcmp(bdy_type, 'Dirichlet'))
    L_h(1,1) = 0; L_h(end,end) = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% OPTION A : CONTINOUS GALERKIN IN TIME %%%%%%%%%%%%%%%%%%%%%%%
%% assembly overall matrix
for ts=1:time_steps
    A_blks{ts} = M_h + ht*L_h;
    B_blks{ts} = -M_h;
    B_blks{ts}(1,1) = 0;
    B_blks{ts}(end,end) = 0;
end

A_blks{1} = eye(pts);

% assemble matrices and vectors
A_large = blkdiag(A_blks{:});
B_large = blkdiag(B_blks{1:end-1});
B_large = [zeros(pts, (time_steps-1)*pts); B_large];
B_large = horzcat(B_large, zeros(time_steps*pts, pts));

comp_mat_CG = A_large + B_large;

rhs = ht* M_h * f;

rhs(1:Nx, 1) = u0;

if(strcmp(bdy_type, 'Dirichlet'))
    rhs(1,:) = left_bdy;
    rhs(end,:) = right_bdy;
end

rhs_CG = reshape(rhs, [(pts*time_steps),1]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% OPTION B : DISCONTINOUS GALERKIN IN TIME %%%%%%%%%%%%%%%%%%%%
%% assembly overall matrix

[comp_mat, L_h, K_t, M_h, M_t, N_t] = set_up_M_st(ht, Nx, Nt, 1, bdy_type);

% create, appropriate f, double entries appropriately for DG

f_double_mat = zeros(2*Nx, Nt-1);
for ts = 1:Nt-1
    f_double_mat(1:Nx,ts) = f(:,ts);
    f_double_mat(Nx+1:end, ts)= f(:,ts+1);
end

rhs_mat = kron(M_t, M_h)*f_double_mat;
% TODO: how to properly deal with initial conditions ... 
rhs_mat(1:Nx, 1) = u0';

if(strcmp(bdy_type, 'Dirichlet'))
    rhs_mat(1, 1:Nt-1) = left_bdy(1:Nt-1);
    rhs_mat(Nx, 1:Nt-1) = right_bdy(1:Nt-1);
    rhs_mat(Nx+1, 1:Nt-1) = left_bdy(2:Nt);
    rhs_mat(end, 1:Nt-1) = right_bdy(2:Nt);
end

rhs_vec = reshape(rhs_mat, [2*Nx*(Nt-1),1]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%spy(tot_mat)
%u_large = tot_mat \ rhs_vec;

u_init_CG = [u0'; zeros(length(rhs_CG)-length(u0), 1)];
u_init_DG = [u0'; zeros(length(rhs_vec)-length(u0), 1)];
max_iter = 10;
smoother = 'GaussSeidel';
eps = 0.25;

u_large_CG = comp_mat_CG \ rhs_CG;
% u_large = comp_mat \ rhs_vec;


%u_large = V_cycle(comp_mat, rhs_vec, u_init, 3, max_iter, smoother);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% SET UP ALGEBRAIC MULTIGRID %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u_large = two_level_AMG(comp_mat, rhs_vec, u_init_DG, smoother, eps, max_iter);

u_mat = zeros(pts, time_steps);
u_mat(:,1) = u0;

u_mat_CG = u_mat;

for ts=2:time_steps
    u_mat_CG(:, ts) = u_large_CG((ts-1)*pts+1:ts*pts);
    u_mat(:, ts) = u_large((2*ts-3)*Nx+1:2*(ts-1)*Nx);
end

figure
mesh(x, t, u_mat')
xlabel('space');
ylabel('time');
zlabel('u');
title('heat equation Neumueller FEM with AMG');

figure
mesh(x, t, u_mat_CG')
xlabel('space');
ylabel('time');
zlabel('u');
title('heat equation FEM w/ CG');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% analysing the coarse and fine grid %%%%%%%%%%%%%%%%%%%

% load('CoarseGridPts.mat');
% C(C == 0) = NaN;
% C(C == 1) = 0;
% 
% F_mat = zeros(Nx, time_steps);
% C_mat = ones(Nx, time_steps);
% for ts=1:time_steps
%     C_mat(:, ts) = C((ts-1)*pts+1:ts*pts);
% end
% 
% [X, T] = meshgrid(x,t);
% 
% figure;
% plot3(X, T, F_mat, 'b.', X, T, C_mat, 'r*');
% title('Fine and Coarse Level Grid Points');
% legend('fine level gridpts', 'coarse level grid pts');
