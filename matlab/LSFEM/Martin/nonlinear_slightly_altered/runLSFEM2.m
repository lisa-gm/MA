% solve nonlinear heat equation on ]0,1[x[0,T] with Neumann b.c.
% nx: spatial grid cells
% nt: time steps
function runLSFEM2(nx,nt)


u = zeros(nx,nt);
v = u;

x = linspace(0,1,nx);
%# u(:,1) = (x<0.75).*(x>0.25); % x>0.5;
%# u(:,1) = sin(2*pi*x);
%u(:,1) = max(0,2*x-1); u = u(:,1)*ones(1,nt);
u(:,1) = max(0,1-2*x); u = u(:,1)*ones(1,nt);

%# u(:,1) = x<0.5;

T = 6;
lsfem(u(:),v(:),nx,nt,1/(nx-1),T/(nt-1),10e-3);
end


  
