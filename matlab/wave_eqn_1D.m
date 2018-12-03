%% Wave Equation %% 
% consider d2u/dt2 = p*d2u/dx2
% rewrite in "laplacian form" and take finite difference scheme
% to build matrix 

%% set up parameters

% for now same number of steps in space and time
N = 10;
delta_x = 1/N;
delta_t = 1/N;

pts = N+1;
x = linspace(0,1, pts);
total_pts = pts * pts;
p = 0.5;

max_it = 10;
u_mat_tot = zeros(max_iter * pts, pts);

max_iter_MG = 30;
levels = 3;

% not including bdy 
%u0 = zeros(1,pts);
u0 = sin(pi*x);
u_init = zeros(total_pts,1);

for i=1:pts
    %u_init((i-1)*pts+1:i*pts) = u0;
end

% for now static
u_bdy_l = 0;
u_bdy_r = 0;

K_h = set_up_K_h_centDiff(N, p, delta_x, delta_t);
%K_h = set_up_K_h_backDiff(N, p, delta_x, delta_t);


% rhs f of changed eqn
f = zeros(total_pts, 1);

% change first and last entry of each block to take bdy conditons
% as well as intial conditions

f(1:pts:end) = u_bdy_l;
f(pts:pts:end) = u_bdy_r;

%% make it a loop to do more time steps
for it=1:max_it
    
f(1:pts) = u0;

% use direct solve for now
u = K_h \ f;
%[u, ~, ~] = W_cycle(K_h, f, u_init, levels, max_iter_MG, 'GaussSeidel');

% reshape u into matrix, each row: one time slice
u_mat = transpose(reshape(u, [pts,pts]));
u_mat_tot((it-1)*pts+1:it*pts,:) = u_mat;

u0 = u_mat(end,:);
end

max_val = max(max(u_mat_tot));
min_val = min(min(u_mat_tot));

%% plot
%figure
    %fh = plot(x, u_mat_tot(1,:));
    %axis([0 1 min_val max_val]);
    %title('1D wave eqn');

for i=2:(max_iter*pts)
    %fh.YData = u_mat_tot(i,:);

    %pause(0.2);
end
    
u_final_t = u_mat_tot(end,:);

size(x)
t= linspace(0,1,size(u_mat_tot,1));
    
figure
mesh(x, t, u_mat_tot)
xlabel('space');
ylabel('time');
zlabel('u');

