clear all;
close all;
beep off;

lw = 1.0;
fs = 12;
ms = 10;

addpath('../src');
addpath('../src/quad_nodes');
addpath('../../mparareal_oo/src');

% Solves y' = lambda*y with Picard iteration
%lambda  = -100;
lambda   = 10*1i;
Tend    = 1.0;

theta_v  = linspace(0,1,500);
nodes_v  = [3 4 5];

gamma    = 1.0;
spec_rad = zeros(3,length(theta_v));

for mm=1:3
    nodes = nodes_v(mm);
    coll    = collocation_sdc(0, Tend, nodes, 'gauss-lobatto');
    for ii=1:length(theta_v)
        theta = theta_v(ii);
        It_mat = coll.getScalarSweepMatrix(lambda, theta, 1.25);
        spec_rad(mm,ii) = max(abs(eig(It_mat)));
    end
end

fig = figure(1);
set(fig, 'Toolbar','none');
set(fig, 'PaperPositionMode', 'auto');
set(fig, 'Units', 'centimeters');
set(fig, 'OuterPosition', [0, 0, 16, 16]);
set(fig, 'Color','white');

plot(theta_v, spec_rad(1,:), 'r', 'linewidth', lw); hold on;
plot(theta_v, spec_rad(2,:), 'b', 'linewidth', lw); 
plot(theta_v, spec_rad(3,:), 'g', 'linewidth', lw);
legend(['M=',num2str(nodes_v(1))], ['M=',num2str(nodes_v(2))], ['M=',num2str(nodes_v(3))]);
title(['\lambda = ', num2str(lambda)]);
xlabel('\theta','fontsize', fs);
ylabel('\rho', 'fontsize', fs);
ylim([0 1.1]);
set(gca, 'fontsize', fs);