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
lambda_v = linspace(-500,-0.001,2000);

Tend    = 1.0;
y0      = 1.0;


nodes        = 3;
q_norm       = zeros(3,1);
q_delta_norm = zeros(3,1);

y_0     = y0*ones(nodes,1);


colls{1}    = collocation_sdc(0, Tend, nodes, 'gauss-legendre');
colls{2}    = collocation_sdc(0, Tend, nodes, 'gauss-lobatto');
colls{3}    = collocation_sdc(0, Tend, nodes, 'gauss-radau-right');

% sol_0   = solution_arraylinear(y0, lambda, 1.0);
% for ii=1:3
%     colls{ii} = colls{ii}.initialize(sol_0);
% end

specrad = zeros(3,length(lambda_v));

for ll=1:length(lambda_v)
    
    lambda = lambda_v(ll);
    
    for ii=1:3
        
        %f       = @(y) lambda*y;
        %y_ex_fh = @(t) y0*exp(lambda*t);
        
        A_sdc = (speye(nodes) - lambda*colls{ii}.Qdelta_mat)\(lambda*(colls{ii}.Qmat - colls{ii}.Qdelta_mat));
        specrad(ii,ll) = max(abs(eig(full(A_sdc))));
       
    end
end

fig = figure(1);
set(fig, 'Toolbar','none');
set(fig, 'PaperPositionMode', 'auto');
set(fig, 'Units', 'centimeters');
set(fig, 'OuterPosition', [0, 0, 16, 16]);
set(fig, 'Color','white');

plot(lambda_v, specrad(1,:), 'b-','linewidth',lw); hold on;
plot(lambda_v, specrad(2,:), 'g-','linewidth',lw);
plot(lambda_v, specrad(3,:), 'r-','linewidth',lw);
legend('Legendre', 'Lobatto', 'Radau','location','southwest');
xlabel('\lambda \Delta t','fontsize',fs);
ylabel('\sigma(M_{sdc})', 'fontsize', fs);
xlim([min(lambda_v) max(lambda_v)]);
ylim([0 1]);
title(['M = ', num2str(nodes)]);
set(gca,'fontsize',fs);

filename = ['sdc_specrad_M', num2str(nodes) '.eps'];
print('-depsc',filename);