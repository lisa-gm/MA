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
lambda  = -5.0;%*1i;
f       = @(y) lambda*y;
Tend    = 1.0;
y0      = 1.0;

y_ex_fh = @(t) y0*exp(lambda*t);

nodes        = 9;
q_norm       = zeros(3,1);
q_delta_norm = zeros(3,1);

y_0     = y0*ones(nodes,1);
sol_0   = solution_arraylinear(y0, lambda, 1.0);

colls{1}    = collocation_sdc(0, Tend, nodes, 'gauss-legendre');
colls{2}    = collocation_sdc(0, Tend, nodes, 'gauss-lobatto');
colls{3}    = collocation_sdc(0, Tend, nodes, 'gauss-radau');
for ii=1:3
   colls{ii} = colls{ii}.initialize(sol_0); 
end

for ii=1:1
    
    A_picard = colls{ii}.Qmat*lambda;
    [EV_pic, D_pic] = eig(A_picard);
    
    A_sdc = (speye(nodes) - lambda*colls{ii}.Qdelta_mat)\(lambda*(colls{ii}.Qmat - colls{ii}.Qdelta_mat));
    
    [EV_sdc, D_sdc] = eig(A_sdc);
    q_norm(ii,1)       = norm(colls{ii}.Qmat);
    q_delta_norm(ii,1) = norm(colls{ii}.Qmat - colls{ii}.Qdelta_mat);
    
    
end

maxmax = max( max(abs(diag(D_pic))), max(abs(diag(D_sdc))) );

% fig = figure(1);
% set(fig, 'Toolbar','none');
% set(fig, 'PaperPositionMode', 'auto');
% set(fig, 'Units', 'centimeters');
% set(fig, 'OuterPosition', [0, 0, 16, 16]);
% set(fig, 'Color','white');
% 
% 
% plot(1:nodes, abs(diag(D_pic)), 'bo-'); hold on;
% plot(1:nodes, abs(diag(D_sdc)), 'ro-');
% legend('Picard','SDC');
% ylim([0 1]);
% xlabel('m');
% ylabel('\lambda_m');
% title(['| \lambda | * dt = ', num2str(abs(lambda)*(Tend-0.0))]);

fig = figure(1);
set(fig, 'Toolbar','none');
set(fig, 'PaperPositionMode', 'auto');
set(fig, 'Units', 'centimeters');
set(fig, 'OuterPosition', [0, 0, 16, 16]);
set(fig, 'Color','white');

subplot(311)
plot(1:nodes, real(EV_pic(:,1:nodes)));
set(gca,'XTick',1:nodes);
ylim([-1 1]);

subplot(312)
plot(1:nodes, imag(EV_pic(:,1:nodes)));
set(gca,'XTick',1:nodes);
ylim([-1 1]);

subplot(313)
plot(1:nodes, abs(diag(D_pic)), 'o-');
ylim([0 maxmax]);

fig = figure(2);
set(fig, 'Toolbar','none');
set(fig, 'PaperPositionMode', 'auto');
set(fig, 'Units', 'centimeters');
set(fig, 'OuterPosition', [0, 0, 16, 16]);
set(fig, 'Color','white');

subplot(311)
plot(1:nodes, real(EV_sdc(:,1:nodes)));
set(gca,'XTick',1:nodes);
ylim([-1 1]);

subplot(312)
plot(1:nodes, imag(EV_sdc(:,1:nodes)));
set(gca,'XTick',1:nodes);
ylim([-1 1]);

subplot(313)
plot(1:nodes, abs(diag(D_sdc)), 'o-');
ylim([0 maxmax]);