
%% space-time solver for heat eq. in 1D
% FEM in space and DG in time


%% time 1D grid (must be equidistant now)
T = 2;
Nt = 40; 
dt = T/(Nt-1);

%% space grid
Nx = 10;
hx = 1/(Nx-1);
x = 0:hx:1;

%u0 = sin(pi*x)';
u0 = 5*ones(Nx,1);


% boundary conditions
bdy_cond = 'Neumann';

%% space operators
Kh = set_up_L_h_FEM(Nx, bdy_cond);
Mh = set_up_M_h_FEM(Nx, bdy_cond);

%% Dirichlet BC
if(strcmp(bdy_cond, 'Dirichlet'))
b_sides = [1 Nx]; 

Kh(b_sides,:) = 0;
Mh(b_sides,:) = 0;

u0(b_sides) = 0;

for i = b_sides
    Kh(i,i) = 1;
    Mh(i,i) = 0;
end

end

%% time operators
syms t

%% basis function on reference element q = 1 
l1(t) = 1-t/dt;
l2(t) = t/dt;
basis = {l1,l2};

K_t = zeros(2);
M_t = zeros(2);
N_tau = zeros(2);

for k = 1:2
    for l = 1:2
        
        bk = basis{k};
        bl = basis{l};
        
        K_t(k,l) = - int(diff(bk,t)*bl,0,dt) + bk(dt)*bl(dt);              
        
        M_t(k,l) = int(bk*bl,0,dt);        
        N_tau(k,l) = bk(0)*bl(dt);       
    end
end

% prova per stabilization
% try! also look at half derivative
% A2 = kron(0.5*(K_t + K_t2),Mh) + kron(M_t,Kh);
A = kron(K_t,Mh) + kron(M_t,Kh);
B = - kron(N_tau,Mh);

Nh = length(A);

tot_size = Nh*(Nt-1);

%% initial condition in time 
f = zeros(1,tot_size);
u0_f = kron([1,0],u0');  % kron([f_space, [1 0]) in C++ il contrario

f(1:2*Nx) = kron([1 0 ; 0 0],Mh) * u0_f'; % questo dovrebbe essere per tutti i blocchi


%% Dirichlet BC
if(strcmp(bdy_cond, 'Dirichlet'))
b_sides = [1 Nx Nx+1 2*Nx]; 

A(b_sides,:) = 0;
B(b_sides,:) = 0;

for i = b_sides
    A(i,i) = 1;
end

end


%% assembling space-time operator
space_time_matrix = zeros(tot_size);

for i = 0:Nt-2
    space_time_matrix(i*Nh+1:(i+1)*Nh,i*Nh+1:(i+1)*Nh) = A;            
     if i > 0
         space_time_matrix(i*Nh+1:(i+1)*Nh,(i-1)*Nh+1:i*Nh) = B;
     end
    
end

% for initial condition in time
% space_time_matrix(1:Nh/2,:) = 0;
% space_time_matrix(1:Nh/2,1:Nh/2) = eye(Nh/2);

space_time_matrix = sparse(space_time_matrix);

%% solve directly
u_dg = space_time_matrix\f';

% transform the solution 
u = [ 0; u0(2:end-1);0]';

for i = Nx+1:2*Nx:length(u_dg) 
    u = [u; transpose(u_dg(i:i+Nx - 1)) ];
end

u = u';
% Other dg node
% for i = 2*Nx+1:2*Nx:length(u_dg)    
%     u = [u; u_dg(i:i+Nx - 1) ];
% end
% u = [u; u_dg(end-Nx+1:end)];

% compute analitical solution
u_true = zeros(size(u));
x = linspace(0,1,Nx);
t = linspace(0,T,Nt);

u0 = 5;

for i = 1:length(x)
   for j = 1:length(t) 
      u_true(i,j) = heat_1D_an( u0,x(i),t(j),Nt);        
   end 
end

u_true(1:Nx,1) = [ 0; u0*ones(Nx-2,1);0]; % initial condition

%% plot
figure
mesh(x, t, u');

error = zeros(1, Nt);

for i=1:Nt
    error(i) = norm(u(:,i) - u_true(:,i));
end

%u(:,1)
%u_true(:,1)
% 
% figure
% plot(0:Nt-1, error);

error;

% plot(u_dg)
% plot(u)
% hold on; plot(u_true)
% grid on 
% legend('DG','true')
% 
% max(abs(u - u_true'))
% max(abs(u(end-Nx/2) - u_true(end-Nx/2)'))
% figure;semilogy(abs(u - u_true')./abs(u_true)','o-')