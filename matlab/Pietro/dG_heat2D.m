%% space-time solver for heat eq. in 2D
% FEM in space and DG in time
clear

addpath('../FEM')
addpath('../meshes')

%% load space grid
mesh = 'structerdSquare4.msh';
mesh_ = read_mesh(mesh);

%% time 1D grid (must be equidistant now)
t0 = 0; T = 1;
N = 2; dt = T/N;
%time = linspace(t0,T,N+1);

p_t = 1; % grado dei polinomi (deve essere 1 per ora)
N_t = p_t + 1; % basis functions per element

%% space operators
[Kh,Mh] = globalStifness(mesh);

%% time operators
% this now works just for p_t = 1
syms t

%% basis function on reference element
l1(t) = 1 - t/dt;
l2(t) = t/dt;
basis = {l1,l2};

K_t = zeros(N_t);
M_t = zeros(N_t);
N_tau = zeros(N_t);

for k = 1:N_t
    for l = 1:N_t
        
        bk = basis{k};
        bl = basis{l};
        
        K_t(k,l) = - int(diff(bk,t)*bl,t0,t0+dt) + bk(t0+dt)*bl(t0+dt);
        
        M_t(k,l) = int(bk*bl,t0,t0+dt);
        
        N_tau(k,l) = bk(t0)*bl(dt);
        
    end
end

A = kron(K_t,Mh) + kron(M_t,Kh);
B = - kron(N_tau,Mh);

Nh = length(A);

% Dirichlet BC TODO(try this before kron!)

b_sides = mesh_.boundary_sides(1,:);

A(b_sides,:) = 0;
A(Nh/N_t+b_sides,:) = 0;

B(b_sides,:) = 0;
B(Nh/N_t+b_sides,:) = 0;


for i = b_sides
    A(i,i) = 1;
    A(Nh/N_t+i,Nh/N_t+i) = 1;
end


%% initial condition in time and BC in space

f = zeros(1,Nh*N);
u0 = 16*mesh_.vertices(1,:).*(1 - mesh_.vertices(1,:)).*(1 - mesh_.vertices(2,:)).*mesh_.vertices(2,:); % initial data

%pdeplot(mesh_.vertices,mesh_.boundary_sides,mesh_.elements,'xydata',u0,'zdata',u0);

% time
f(1:Nh/2) = u0;

%% assembling space-time operator

space_time_matrix = zeros(Nh*N);

for i = 0:N-1
    
    space_time_matrix(i*Nh+1:(i+1)*Nh,i*Nh+1:(i+1)*Nh) = A;
    
    if i > 0
        space_time_matrix(i*Nh+1:(i+1)*Nh,(i-1)*Nh+1:i*Nh) = B;
    end
    
end

space_time_matrix = sparse(space_time_matrix);

cond(space_time_matrix)


% for initial condition in time
% space_time_matrix(1:Nh/2,:) = 0;
% space_time_matrix(1:Nh/2,1:Nh/2) = eye(Nh/2);

%% direct solver
% tic
% u = space_time_matrix\f';
% toc
% 
% plot(u)

% u = zeros(length(f), 1);
% 
% iters = 100;
% error = zeros(iters,1);
% block_size = 3;
% tic;
% %line-search
% for it = 1:iters
%     u =  block_GS(u, -1, f', space_time_matrix, 1, block_size,0);
%     error(it) = norm(space_time_matrix*u - f');
%     if error(it) < 1e-14
%         break;
%     end
% end
% toc;
% 
% %standard
% error2 = zeros(iters,1);
% u = zeros(length(f), 1);
% 
% tic;
% for it = 1:iters
%    u =  block_GS(u, 1, f', space_time_matrix, 1, block_size,0);
%    error2(it) = norm(space_time_matrix*u - f');
%    
%      if error2(it) < 1e-14
%         break;
%     end
% end
% toc;
% 
% close all;
% semilogy(error); 
% hold on;
% semilogy(error2, 'r'); 

% %% movie 
% % 
% % %Preallocate movie structure.
% M(1:N*N_t) = struct('cdata', [],'colormap', []);
% 
% for n = 1:N*N_t
%     pdeplot(mesh_.vertices,mesh_.boundary_sides,mesh_.elements,'xydata',u(1+(n-1)*Nh/2:n*Nh/2),'zdata',u(1+(n-1)*Nh/2:n*Nh/2));    
%     xlim([0 1]);
%     ylim([0 1]);
%     zlim([0 1]);
%     drawnow
%     M(n) = getframe;
% end
% movie(M);


