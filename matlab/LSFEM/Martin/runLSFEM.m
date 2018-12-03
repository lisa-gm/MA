 % solve standard heat equation on ]0,1[^2 with Neumann b.c.
% nx: spatial grid cells
% nt: time steps
function runLSFEM(nx,nt)


u = zeros(nx,nt);
v = u;

x = linspace(0,1,nx);
u(:,1) = (x<0.75).*(x>0.25); % x>0.5;
u(:,1) = sin(2*pi*x);


lsfem(u(:),v(:),nx,nt,1/(nx-1),1/(nt-1));
