%% advection equation (*) du/dt + a * du/dx = 0
% solve numerically after rewriting derivatives
% in direction of characteristics and orthogonal part
% characteristics: x = at + x0
% the vectors are v_n = 1/sqrt(1+a^2)*(-1 a)
% v_t = 1/sqrt(1+1/a^2)*(1 1/a)

% we obtain the following equation
% (*) = (a*sqrt(1+1/a^2)/(a^2+1) + a^3*sqrt(1+1/a^2)/(a^2+1))*du/dvt = 0
% du/dvn part canceled out
% hence we see that u is constant along vt

%% theoretical grid set up
% did we assume an equidistant grid?
% have to be careful to match up the right grid points with each other

% how do we make sure we have grid points exactly pointing in the direction
% of vt?
% can't rotate grid, because then we would also rotate boundary conditions

% need to get the ratio of hx and ht right
% suppose we fix hx,  how do we then have to choose ht?
% we have t = 1/a * x - x0/a
% hence if hx = 1/(n-1) then ht = 1/a*hx and then we can use:
% (u_ij - u_i-1j-1) / tau = 0 

%% setting up the parameters
clear all;
close all;

a = 2;

n = 30;
x = linspace(0,1, n);
hx = x(2);
%ht = 1/a * hx;
n_t = 80;
ht = 1/(n_t -1);

CFL = a*ht/hx;
fprintf('CFL = %d\n', CFL);

T = 1;
t = 0:ht:T;
time_steps = length(t);

% BC and IC
u0 = sin(4*pi*x);
%u0 = ones(n,1);
u0(ceil(n/4):end) = 0;

%figure
%plot(x, u0);

bdy_left = zeros(time_steps, 1);
bdy_right = zeros(time_steps, 1);


%% set up matrix and rhs 

pts = time_steps*n;
temp = ones(pts,1);
% 1/ht or do i use the distance between the actual points? 
dist = (a^2 + 1) / (abs(a)*(n-1));
factor = 1/dist*(a*sqrt(1+1/a^2)/(a^2+1) + a^3*sqrt(1+1/a^2)/(a^2+1));
M = factor*spdiags([temp, -temp], [0, -(n+1)], pts, pts);
M(1:n, 1:n) = eye(n,n);
M(1:n:end, :) = 0;
M(n:n:end, :) = 0;
M(1:n:end, 1:n:end) = eye(time_steps, time_steps);
M(n:n:end, n:n:end) = eye(time_steps, time_steps);

%spy(M)

rhs = zeros(pts, 1);
rhs(1:n) = u0;
rhs(1:n:end) = bdy_left;
rhs(n:n:end) = bdy_right;

%% solve

u = M \ rhs;

%% plot results

u_mat = zeros(time_steps, n);
for ts=1:time_steps
    u_mat(ts, :) = u((ts-1)*n+1:ts*n);
end
% 
figure
mesh(x, t, u_mat)
xlabel('space');
ylabel('time');
zlabel('u');
title('advection equation computed along characteristics');

% how do I test now if this works better than a solver that is aligned with
% the coordinate axis? 

% compute actual solution
% compute solution with axis aligned solver
% ie (u_ij - u_ij-1) / ht + a (u_ij - u_i-1j) / hx = 0

%% set up matrix 

n_st = n;
hx_st = 1/(n_st-1);
x_st = linspace(0,1, n_st);

T = 1;
t_st = linspace(0,T, time_steps);
ht_st = t_st(2);
ts_st = length(t_st);

bdy_left_st = zeros(ts_st,1);
bdy_right_st = zeros(ts_st,1);

u0_st = sin(4*pi*x_st);

%u0_st = ones(n_st, 1);
u0_st(ceil(n_st/4):end) = 0;

pts = ts_st*n_st;

rhs_st = zeros(pts,1);
rhs_st(1:n_st) = u0_st;
rhs_st(1:n_st:end) = bdy_left_st;
rhs_st(n:n_st:end) = bdy_right_st;


temp = ones(pts,1);
M_st = spdiags([1/ht_st*temp, (-1/ht+a/hx_st)*temp, -a/hx_st*temp], [0, -n_st, -(n_st+1)], pts, pts);
M_st(1:n_st, 1:n_st) = eye(n_st,n_st);
M_st(1:n_st:end, :) = 0;
M_st(n_st:n_st:end, :) = 0;
M_st(1:n_st:end, 1:n_st:end) = eye(ts_st, ts_st);
M_st(n_st:n_st:end, n_st:n_st:end) = eye(ts_st, ts_st);

%spy(M_st)

%% solve

u_st = M_st \ rhs_st;

u_mat_st = zeros(ts_st, n_st);
for ts=1:ts_st
    u_mat_st(ts, :) = u_st((ts-1)*n_st+1:ts*n_st);
end

figure
mesh(x_st, t_st, u_mat_st)
xlabel('space');
ylabel('time');
zlabel('u');
title('advection equation computed along axes');


%% compute analytical solution 
% by shifting it across the board

u_an = zeros(time_steps, 2*n);
for i=1:min(n, time_steps)
    u_an(i, i:i+n-1) = u0';
end

% cut off extra bits
u_an = u_an(:, 1:n);

% size(x)
% size(t)
% size(u_an)

figure
mesh(x, t, u_an)
xlabel('space');
ylabel('time');
zlabel('u');
title('analytical solution');

error_charact_sol = norm(u_an(:,1:end-1) - u_mat(:,1:end-1));
error_normal_sol = norm(u_an(:,1:end-1) - u_mat_st(:,1:end-1)); 

fprintf('error_charact_sol = %d \n', error_charact_sol);
fprintf('error_normal_sol = %d \n', error_normal_sol);

