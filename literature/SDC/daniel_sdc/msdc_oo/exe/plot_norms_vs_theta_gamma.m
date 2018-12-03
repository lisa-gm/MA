clear all;
close all;
clear classes;
beep off;
addpath('../src');
addpath('../src/quad_nodes');
addpath('../../mparareal_oo/src');
ms = 10;
fs = 20;
lw = 1.5;

M = 3;
gamma_v = linspace(0.6, 2.0, 200);
theta_v = linspace(0.0, 2.0, 200);
lambda = 0 + 10*1i;

% type = 'gauss-radau-right';
type = 'gauss-legendre';

E_sdc = zeros(length(theta_v), length(gamma_v));


for mm=1:length(theta_v)
    coll = collocation_sdc(0, 1, M, type);
    for ll=1:length(gamma_v)
        Mit          = coll.getScalarSweepMatrix( lambda, theta_v(mm), gamma_v(ll) );
        E_sdc(mm,ll)   = max(abs(eig(full(Mit))));
    end
    
end

color{1} = [0.65 0.65 0.65];
color{2} = 'k';


fig = figure(1);
set(fig, 'Toolbar','none');
set(fig, 'PaperPositionMode', 'auto');
set(fig, 'Units', 'centimeters');
set(fig, 'OuterPosition', [0, 0, 16, 16]);
set(fig, 'Color','white');
contour(gamma_v, theta_v, E_sdc, 0.999, 'linewidth', lw+0.5, 'color', color{1}); hold on;
[C1,h1] = contour(gamma_v, theta_v, E_sdc, 0.1:0.1:0.9, 'linewidth', lw, 'color', color{2});
set(h1, 'showtext','on');
set(gca,'fontsize',fs);
xlim([min(gamma_v) max(gamma_v)]);
ylim([min(theta_v) max(theta_v)]);
xlabel('\gamma', 'fontsize', fs);
ylabel('\theta', 'fontsize', fs);

if real(lambda)==0
    lambda_text = 'imag';
else
    lambda_text = 'real';
end
% 
filename = ['specrad_', lambda_text '_', num2str(abs(lambda)) '_M', num2str(M) '.eps'];
print('-depsc',filename);