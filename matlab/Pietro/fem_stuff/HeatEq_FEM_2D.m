clear all

[A,Mass] = globalStifness('Lshape.msh');

image(Mass)

T_final = 0.4;
dt = 0.1;
t_steps = 1 + T_final/dt;

mesh = read_mesh('Lshape.msh');
u_0 = mesh.vertices(1,:).*mesh.vertices(2,:); % initial data

% prealloc
u = zeros(T_final/dt,length(u_0));
u(1,:) = u_0;

% no BC imposed <==> Neumann BC (a.k.a. nataural BC)

for t = 1:t_steps
    u(t+1,:)=(Mass+dt*A)\Mass*u(t,:)'; 
end

%--------figure editing---------%
z_max = max(max(u_0));

subplot(2,2,1);pdemesh(mesh.vertices,mesh.boundary_sides,mesh.elements,u(1,:));zlim([0,z_max]);grid on;
title('initial time');
subplot(2,2,2);pdemesh(mesh.vertices,mesh.boundary_sides,mesh.elements,u(floor(t_steps/2),:));zlim([0,z_max]);grid on;
title('middle time');
subplot(2,2,3);pdemesh(mesh.vertices,mesh.boundary_sides,mesh.elements,u(t_steps,:));zlim([0,z_max]);grid on;
title('final time');

% invariance of the integral
S=zeros(1,t_steps);
for n=1:t_steps
    S(n)=sum(Mass*u(n,:)');
end
subplot(2,2,4);plot(S,'.-');
title 'integral of u(x) in each time step'
    
    
%----------movie editing-----------%    
% % Preallocate movie structure.
%  M(1:T_final/tau) = struct('cdata', [],'colormap', []);
% 
% for n=1:T_final/tau
% pdeplot(mesh.vertices,mesh.boundary_sides,mesh.elements,'xydata',u(n,:),'zdata',u(n,:));
% axis tight;
% axis equal;
% M(n)=getframe;
% end
% movie(M);

