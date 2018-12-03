clear all;
close all;
clear classes;
beep off;
addpath('../src');
addpath('../src/quad_nodes');
addpath('../../mparareal_oo/src');
ms = 10;
fs = 12;
lw = 1.5;

Mvec = 2:30;
theta = 1.0;
gamma = 1.0;
lambda = -100;
type = {'gauss-legendre', 'gauss-radau-right', 'gauss-chebyshev'};

E_sdc = zeros(length(type), length(Mvec));
est_sdc = zeros(length(type), length(Mvec));
sr_sdc = zeros(length(type), length(Mvec));

for tt=1:length(type)
    
    for ll=1:length(Mvec)
        
        coll           = collocation_sdc(0, 1, Mvec(ll), type{tt});
        %Mit            = coll.getScalarSweepMatrix(lambda, theta, gamma);
        Mit            = coll.Qmat - coll.Qdelta_mat;
        Mprecond       = (speye(Mvec(ll)) - lambda*coll.Qdelta_mat)\speye(Mvec(ll));
        est_sdc(tt,ll) = coll.lebesgue;
%         E_sdc(tt,ll)   = norm(Mprecond, inf);
        E_sdc(tt,ll)   = norm(coll.Qmat, inf);
        
    end
    
end

fig = figure(1);
set(fig, 'Toolbar','none');
set(fig, 'PaperPositionMode', 'auto');
set(fig, 'Units', 'centimeters');
set(fig, 'OuterPosition', [0, 0, 16, 16]);
set(fig, 'Color','white');

plot(Mvec, E_sdc(1,:), 'bo', 'linewidth', lw); hold on;
plot(Mvec, E_sdc(2,:), 'rd', 'linewidth', lw);
plot(Mvec, E_sdc(3,:), 'gs', 'linewidth', lw);
legend(type{1}, type{2}, type{3}, 'location','northwest');
plot(Mvec, est_sdc(1,:), 'b-', 'linewidth', lw);
plot(Mvec, est_sdc(2,:), 'r-', 'linewidth', lw);
plot(Mvec, est_sdc(3,:), 'g-', 'linewidth', lw);
plot(Mvec, 1.0+0*E_sdc(1,:), 'k-');
xlabel('M','fontsize',fs);
ylabel('|| Q ||','fontsize',fs);
print('-depsc','norm_plot.eps');

fig = figure(2);
set(fig, 'Toolbar','none');
set(fig, 'PaperPositionMode', 'auto');
set(fig, 'Units', 'centimeters');
set(fig, 'OuterPosition', [0, 0, 16, 16]);
set(fig, 'Color','white');

plot(Mvec, E_sdc(1,:), 'bo', 'linewidth', lw); hold on;
plot(Mvec, E_sdc(2,:), 'rd', 'linewidth', lw);
plot(Mvec, E_sdc(3,:), 'gs', 'linewidth', lw);
legend(type{1}, type{2}, type{3}, 'location','northwest');
plot(Mvec, 1.0+0*E_sdc(1,:), 'k-');
xlabel('M','fontsize',fs);
ylabel('|| Q ||','fontsize',fs);
print('-depsc','norm_plot_2.eps');

