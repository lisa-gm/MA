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

Mvec = [3 3 3];
lambda_v = linspace(-100, -1e-2, 1e4);
theta = [1.0 1.0 1.0];
gamma = [1.0 1.15 1.3];
%type = 'gauss-chebyshev';
% type = {'gauss-legendre', 'gauss-radau', 'gauss-chebyshev'};
% type = {'gauss-chebyshev', 'gauss-chebyshev', 'gauss-chebyshev'};
% type = 'gauss-radau-right';
type = 'gauss-legendre';

E_sdc = zeros(length(Mvec), length(lambda_v));
est_sdc = zeros(length(Mvec), length(lambda_v));
sr_sdc = zeros(length(Mvec), length(lambda_v));

for mm=1:length(Mvec)
    
    coll = collocation_sdc(0, 1, Mvec(1), type);
    
    for ll=1:length(lambda_v)
        lambda = 1i*lambda_v(ll);
        
        Mit          = coll.getScalarSweepMatrix( lambda, theta(mm), gamma(mm));
        %Mright   = lambda*( coll.Qmat - coll.Qdelta_mat);
        % Mprecond = (speye(Mvec(mm)) - lambda*coll.Qdelta_mat)\speye(Mvec(mm));
        %cc = 1 - lambda*coll.delta_m;
        %precond = sum(2./abs(cc));
        est_sdc(mm,ll) = abs(lambda)*(2 + coll.lebesgue );
        E_sdc(mm,ll)   = max(abs(eig(full(Mit))));
        %E_sdc(mm,ll)   = norm(Mit, inf);
        sr_sdc(mm,ll) = max(abs(eig(full(Mit))));
    end
    
end

fig = figure(1);
set(fig, 'Toolbar','none');
set(fig, 'PaperPositionMode', 'auto');
set(fig, 'Units', 'centimeters');
set(fig, 'OuterPosition', [0, 0, 16, 16]);
set(fig, 'Color','white');

semilogx(lambda_v, E_sdc(1,:), 'b-', 'linewidth', lw); hold on;
semilogx(lambda_v, E_sdc(2,:), 'r-', 'linewidth', lw);
semilogx(lambda_v, E_sdc(3,:), 'g-', 'linewidth', lw);
% legend(['M=',num2str(Mvec(1))], ['M=',num2str(Mvec(2))], ['M=',num2str(Mvec(3))], 'location', 'northwest'); legend boxoff;
legend(['\gamma=',num2str(gamma(1))], ['\gamma=',num2str(gamma(2))], ['\gamma=',num2str(gamma(3))], 'location', 'northwest'); legend boxoff;
% semilogy(lambda_v, est_sdc(1,:), 'b--'); hold on;
% semilogy(lambda_v, est_sdc(2,:), 'r--');
% semilogy(lambda_v, est_sdc(3,:), 'g--');
ylim([0 0.8]);
xlim([min(lambda_v) max(lambda_v)]);

% fig = figure(2);
% set(fig, 'Toolbar','none');
% set(fig, 'PaperPositionMode', 'auto');
% set(fig, 'Units', 'centimeters');
% set(fig, 'OuterPosition', [0, 0, 16, 16]);
% set(fig, 'Color','white');
% 
% plot(lambda_v, sr_sdc(1,:), 'b-', 'linewidth', lw); hold on;
% plot(lambda_v, sr_sdc(2,:), 'r-', 'linewidth', lw);
% plot(lambda_v, sr_sdc(3,:), 'g-', 'linewidth', lw);
% legend(['M=',num2str(Mvec(1))], ['M=',num2str(Mvec(2))], ['M=',num2str(Mvec(3))], 'location', 'southwest'); legend boxoff;