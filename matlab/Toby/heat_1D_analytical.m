%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% checking derivatives %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%clear all;
%close all;

S = 1;
Nx = 40;
hx = 1/(Nx-1);
x = linspace(0,S, Nx);

T = 3;
Nt = 30;
ht = 1/(Nt-1);
t = linspace(0,T, Nt);

u_mat = zeros(Nx, Nt);

u0=3;
n_max = 20;


for j=1:Nt
    for i=1:Nx
        u_mat(i,j) = heat_1D_analyt(u0, x(i), t(j), n_max);
    end
end

figure
mesh(x, t, u_mat');


function u = heat_1D_analyt(u0,x,t,n_max)

% u0 : initial (constant) value
% n : bigger n better solution 

u = 0;

for n = 1:n_max 
    N = 2*n-1;
    u = u + sin(N*pi*x)*exp(-N^2*pi^2*t)/N;    
end

u =  (4*u0/pi) * u;

end