%% implicit wave equation, small matrices

% of the form d^2u/dt^2-c^2d^2u/dx^2 = f
clear all;
close all;

c2 = 4;
c = sqrt(c2);

int_time = 10;
int_space = 1;

Nx = 50;
Nt = 1000;

hx = int_space/(Nx-1);
ht = int_time/(Nt-1);

x = linspace(0,int_space,Nx);
t = linspace(0,int_time,Nt);


% add in f, for now time independent
%f = 0.02*cos(pi*x);
%f = zeros(Nx+1,1);

u0 = sin(pi*x);
%u0(10:end) = 0;

u = zeros(Nt, Nx);
u(1,:) = u0;

u1 = zeros(Nx,1);

u(1,:) = u0';
u(2,:) = u0';

r = -c2*ht^2/hx^2

temp = ones(Nx,1);
L_h = spdiags([r*temp, -2*r*temp, r*temp], -1:1, Nx, Nx);

% zero boundary conditions PUT BACK IN
L_h(1,:) = [0, zeros(1,Nx-1)];
L_h(end,:) = [zeros(1,Nx-1), 0];

I = eye(Nx, Nx);

for ts = 3:Nt
   % if boundary not zero, need to make changs here...
   u(ts,:) = (I + L_h) \ (2*u(ts-1,:)' - u(ts-2, :)');
  %u(ts, :) = L_h*u(ts-1,:)'-u(ts-2,:)';
end

(I + L_h)

figure
mesh(x, t, u)
xlabel('space');
ylabel('time');
zlabel('u');
title('wave equation implicit scheme');

%u(end, :)
%-sin(2*pi*x)
%figure 
for ts=1:Nt+1
%plot(x, u(ts,:));
%axis([0 1 -1 1]);

%pause(0.15);
end