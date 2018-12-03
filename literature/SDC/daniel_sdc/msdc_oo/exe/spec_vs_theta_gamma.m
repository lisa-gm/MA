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
lambda_v= linspace(10, 1, 50)*(-1 + 1i);
Tend    = 1.0;

theta_v  = linspace(0,1,100);
gamma_v  = linspace(0.0,2.0,200);
nodes    = 9;

spec_rad = zeros(length(gamma_v),length(theta_v));
coll    = collocation_sdc(0, Tend, nodes, 'gauss-legendre');

for ll=1:length(lambda_v)
    lambda = lambda_v(ll);
    for jj=1:length(gamma_v)
        gamma = gamma_v(jj);
        for ii=1:length(theta_v)
            theta = theta_v(ii);
            It_mat = ( speye(nodes) - theta*lambda*coll.Qdelta_mat )\( (1 - gamma)*speye(nodes) + lambda*(coll.Qmat - theta*coll.Qdelta_mat) );
            spec_rad(jj,ii) = max(abs(eig(It_mat)));
        end
    end
    
    if ll==1
        s_max = max(max(spec_rad));
        s_min = min(min(spec_rad));
    end
    
    fig = figure(1); clf;
    set(fig, 'Toolbar','none');
    set(fig, 'PaperPositionMode', 'auto');
    set(fig, 'Units', 'centimeters');
    set(fig, 'OuterPosition', [0, 0, 16, 16]);
    set(fig, 'Color','white');
    
    surface(theta_v, gamma_v, spec_rad, 'LineStyle', 'none', 'EdgeColor', 'none'); hold on;
    shading interp;
    colorbar;
    [xmin, imin] = min(spec_rad);
    [ymin, jmin] = min(xmin);
    plot3(theta_v(jmin),gamma_v(imin(jmin)),2*ymin,'ro','MarkerSize',10, 'markerfacecolor','r');
    xlabel('\theta','fontsize',fs);
    ylabel('\gamma','fontsize',fs);
    caxis([s_min s_max]);
    title(['\lambda = ', num2str(lambda, '%3.2f') ' -- M = ', num2str(nodes)]);
    drawnow;
    
    if ll<10
        filename = ['theta_gamma_00',num2str(ll) '.png'];
    elseif ll<100
        filename = ['theta_gamma_0',num2str(ll) '.png'];
    else
        filename = ['theta_gamma_',num2str(ll) '.png'];
    end
    print('-dpng','-r300',filename);
end