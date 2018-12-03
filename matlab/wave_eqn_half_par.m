%% try more explicit scheme
clear all;
%close all;

% what I want to solve here, u_ij where i is space, j time
% 1/dt^2*(u_ij - 2u_ij-1 + u_ij-2) - c^2/dx^2 (u_i+1j-1 -2uij-1 + ui-1j-1)
% = f_ij

%% set up parameters
c2 = 4;
c = sqrt(c2);

int_time = 1;
int_space = 1;

Nx = 50;
% breaks for Nt = 500, k = c^2*dt^2/dx^2 = 4*1/99^2*49^2
% only stable for values <1, but this is the case ..?!
% .... why?!
Nt = 100;

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
max(u0)
%u0(10:end) = 0;

% assemble rhs
rhs = zeros(Nx*Nt,1);

% put in initial conditions
rhs(1:Nx) = u0;
rhs(Nx+1:2*Nx) = u0;

% put in boundary conditions
% still have to do scaling?!
rhs(1:Nx:end) = bdy_l;
rhs(Nx:Nx:end) = bdy_r;

pts = Nx;
total_pts = Nt*Nx;   

%% set up matrix
% build B as a diagonal matrix first then add zeros at the top and right side
% to make it lower diagonal, see structure through spy
    temp = ones(total_pts,1);

    B = spdiags([-c2/hx^2*temp, -2*(1/ht^2-c2/hx^2)*temp, -c2/hx^2*temp], -1:1, pts, pts);
    B_list = {};
    
    for i=1:Nt-1
        B_list{i} = B;
    end
    
    K_h_blk = blkdiag(B_list{1:end});
   
    % this is actually for lower off diagonal blocks
    % have to add on to make the matrix large enough
    % and then put ones on diagonal
    
    K_h_off_blk = [zeros(pts, total_pts-pts); K_h_blk];
    
    % now glue zeros vector on the side
    K_h_off_blk = horzcat(K_h_off_blk, zeros(total_pts, pts));
    
    % get the diagonal and lower lower diagonal entries
    K_h_ones =  spdiags([1/ht^2*temp, 1/ht^2*temp], [0,-(2*pts)], total_pts, total_pts);
    K_h = (K_h_off_blk + K_h_ones);
        
    % the first 2 blocks cannot change, we have the initial conditions,
    % hence replace first row block by identity and then zeros
    
    K_h(1:2*pts,1:2*pts) = eye(2*pts);
    K_h(1:2*pts,2*pts+1:end) = zeros(2*pts, total_pts - 2*pts);
    
    % we also have fixed boundary conditions for first and last entry
    % of every block in vector, make rows with entry 1
    
    K_h(1:pts:end, :) = 0;
    K_h(pts:pts:end, :) = 0; 
    
    K_h(1:pts:end, 1:pts:end) = eye(Nt);
    K_h(pts:pts:end, pts:pts:end) = eye(Nt); 
    
    %figure
    %spy(K_h)
    
%% solve 

u = K_h \ rhs;
    
%% get result into matrix format

u_mat = reshape(u, [Nx, Nt]);
u_mat = u_mat';
% norm(u_mat(1, :) - sin(2*pi*x))

%% plot result
figure
mesh(x, t, u_mat);
xlabel('space');
ylabel('time');
zlabel('u');
title('wave eqn half exp');
% 
% u_analyt = zeros(Nt, Nx);
% for ts=1:Nt
%        for sp=1:Nx
%     u_analyt(ts, sp) = sin(pi*x(sp))*cos(2*pi*t(ts));
%        end
% end

max(max(u_mat))
min(min(u_mat))

% figure
% mesh(x, t, u_analyt);
% xlabel('space');
% ylabel('time');
% zlabel('u');
% title('analytical solution');
% 
% u0_a = u_analyt(1,:);
% u1_a = u_analyt(2,:);
% u2_a = u_analyt(3,:);
