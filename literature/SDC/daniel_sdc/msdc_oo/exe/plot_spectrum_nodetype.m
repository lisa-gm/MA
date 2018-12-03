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
lambda_v = [-100 -1];

Tend    = 1.0;
y0      = 1.0;


nodes        = 9;
q_norm       = zeros(3,1);
q_delta_norm = zeros(3,1);

y_0     = y0*ones(nodes,1);


colls{1}    = collocation_sdc(0, Tend, nodes, 'gauss-legendre');
colls{2}    = collocation_sdc(0, Tend, nodes, 'gauss-lobatto');
colls{3}    = collocation_sdc(0, Tend, nodes, 'gauss-radau');

specrad = zeros(2,3,nodes);

for ll=1:length(lambda_v)
    
    lambda = lambda_v(ll);
    
    for ii=1:3
        
        %f       = @(y) lambda*y;
        %y_ex_fh = @(t) y0*exp(lambda*t);
        
        A_sdc = (speye(nodes) - lambda*colls{ii}.Qdelta_mat)\(lambda*(colls{ii}.Qmat - colls{ii}.Qdelta_mat));
        specrad(ll,ii,:) = eig(full(A_sdc));

    end
end

for ll=1:length(lambda_v)
    
    fig = figure(ll);
    set(fig, 'Toolbar','none');
    set(fig, 'PaperPositionMode', 'auto');
    set(fig, 'Units', 'centimeters');
    set(fig, 'OuterPosition', [0, 0, 16, 16]);
    set(fig, 'Color','white');
    
    plot( reshape(real(specrad(ll,1,:)), nodes, 1), reshape(imag(specrad(ll,1,:)), nodes, 1), 'bo', 'markerfacecolor', 'b', 'markersize', ms); hold on;
    plot( reshape(real(specrad(ll,2,:)), nodes, 1), reshape(imag(specrad(ll,2,:)), nodes, 1), 'gd', 'markerfacecolor', 'g', 'markersize', ms);
    plot( reshape(real(specrad(ll,3,:)), nodes, 1), reshape(imag(specrad(ll,3,:)), nodes, 1), 'rs', 'markerfacecolor', 'r', 'markersize', ms);
    legend('Legendre','Lobatto','Radau','location','northwest');

    xlabel('Real','fontsize',fs);
    ylabel('Imag', 'fontsize', fs);
    title(['\lambda = ', num2str(lambda_v(ll)) ', M = ', num2str(nodes)]);
    
    filename = ['sdc_spectrum_lambda', num2str(abs(lambda_v(ll))) '_M', num2str(nodes) '.eps'];
    print('-depsc',filename);
    
end