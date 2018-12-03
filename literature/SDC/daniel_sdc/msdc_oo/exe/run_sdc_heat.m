clear all;
close all;
beep off;
clear classes;

lw = 1.4;
fs = 12;
ms = 10;

addpath('../src');
addpath('../src/quad_nodes');
addpath('../../mparareal_oo/src');
%addpath('../../matrix_mg');
addpath('~/Programs/Matlab/ExampleCodes/LinearParareal');

Nx = 201;
nu = 10.0;
%[A, mesh] = getLapMat(0, 1, Nx);
[A, mesh] = getLaplace(0, 1, Nx, 4);
A         = nu*A;

y0        = sin(2*pi*mesh).';
y_ex_fh   = @(t) exp(-4*nu*pi^2*t)*y0;

Tend   = 1/50;
Nsteps = 1;
Tmesh  = linspace(0,Tend,Nsteps+1);
nodes  = 7;
Nit    = 25;

gamma = 1.0;

method = {'full', 'full'};
Nit_smoother = 1;

sol_0   = solution_arraylinear(y0,            A, speye(Nx));
%sol_0   = solution_arraylinear_iter(y0,       A, speye(Nx), 0, 1, nu, method, 1);
sol_end = solution_arraylinear(y_ex_fh(Tend), A, speye(Nx));

eigv = zeros(1,3);

colls{1}    = collocation_sdc(0, Tend, nodes, 'gauss-legendre');
colls{2}    = collocation_sdc(0, Tend, nodes, 'gauss-radau-right');
%colls{3}    = collocation_sdc(0, Tend, nodes, 'gauss-radau');
%colls{4}    = collocation_sdc(0, Tend, nodes, 'equi-noleft');

nr_types = length(colls);

for ii=1:nr_types
   colls{ii} = colls{ii}.initialize(sol_0);
end


% For each set of nodes, generate cell array that contains the exact
% solution at the nodes
sol_ex      = cell(3,nodes);
for ii=1:nr_types
    for nn=1:nodes
        sol_ex{ii,nn} = solution_arraylinear(y_ex_fh(colls{1,ii}.nodes(nn)), A, speye(Nx));
    end
end

% For each set of nodes, compute also the correct collocation solution by
% solving U - Q*lambda*U = U0
sol_coll_ex = cell(nr_types,nodes);
for ii=1:nr_types
    
    u_coll = (speye(nodes*Nx) - kron(colls{1,ii}.Qmat,A))\(kron(ones(nodes,1),y0));
    u_coll = reshape(u_coll, Nx, nodes);
    for nn=1:nodes
        
        sol_coll_ex{ii,nn} = solution_arraylinear( u_coll(:,nn), A, speye(Nx));
    end
end

error    = zeros(nr_types,Nit);
err_end  = zeros(nr_types,Nit);
residual = zeros(nr_types,Nit);

% Now perfom Nit iterations
for nn=1:Nit
    for ii=1:nr_types
        
        % Do Picard iteration
        colls{ii} = colls{ii}.applySweep(1.0, gamma);
        
        % Fetch approximation at Tend
        uend            = colls{ii}.get_end_value;
        
        % Compute difference to analytical solution and compute error
        diff            = uend.axpy(-1.0, sol_end);
        err_end(ii,nn)  = diff.getnorm;
        
        % Fetch residual
        residual(ii,nn) = colls{ii}.get_residual;
        
        % compute defect on all nodes to analytical solution
        def = zeros(nodes,1);
        for mm=1:nodes
            diff      = sol_coll_ex{ii,mm}.axpy(-1.0, colls{ii}.node_solutions{1,mm});
            def(mm,1) = diff.getnorm;
        end
        
        error(ii,nn) = norm(def, inf);
    end
    
end

res_tol = 1e-6;
it_res  = find(residual(1,:)<res_tol, 1, 'first');
fprintf('Iterations until residual tolerance: %i \n', it_res);
fprintf('Resulting error: %5.2e \n', error(1,it_res));

fig = figure(1);
set(fig, 'Toolbar','none');
set(fig, 'PaperPositionMode', 'auto');
set(fig, 'Units', 'centimeters');
set(fig, 'OuterPosition', [0, 0, 16, 16]);
set(fig, 'Color','white');
semilogy(1:Nit, reshape(error(1,:), 1, Nit), 'ro:', 'linewidth', lw, 'markersize', ms); hold on;
semilogy(1:Nit, reshape(error(2,:), 1, Nit), 'bo:', 'linewidth', lw, 'markersize', ms);
% semilogy(1:Nit, reshape(error(3,:), 1, Nit), 'go:', 'linewidth', lw, 'markersize', ms);
legend('Gauss-Legendre','Gauss-Radau', 'location','northeast'); legend boxoff;
% if nr_types==4
%    semilogy(1:Nit, reshape(error(4,:), 1, Nit), 'kx:', 'linewidth', lw, 'markersize', ms);
%    legend('Gauss-Legendre','Gauss-Lobatto','Gauss-Radau','Equi-noleft','location','northeast'); legend boxoff;
% else
%     legend('Gauss-Legendre','Gauss-Lobatto','Gauss-Radau','location','northeast'); legend boxoff;
% end
ylim([1e-16 1e2]);
xlim([1 Nit]);
ylabel('Iteration error','fontsize',fs);
xlabel('Iteration','fontsize',fs);
title(['N_x = ', num2str(Nx), ', Smoother : ', method{1} '(',num2str(Nit_smoother) '),  M=', num2str(nodes)]);
set(gca,'fontsize',fs);

% fig = figure(1);
% set(fig, 'Toolbar','none');
% set(fig, 'PaperPositionMode', 'auto');
% set(fig, 'Units', 'centimeters');
% set(fig, 'OuterPosition', [0, 0, 32, 24]);
% set(fig, 'Color','white');
% 
% subplot(131)
% semilogy(1:Nit, reshape(error(1,:), 1, Nit), 'ro:', 'linewidth', lw, 'markersize', ms); hold on;
% semilogy(1:Nit, reshape(error(2,:), 1, Nit), 'bo:', 'linewidth', lw, 'markersize', ms);
% semilogy(1:Nit, reshape(error(3,:), 1, Nit), 'go:', 'linewidth', lw, 'markersize', ms);
% if nr_types==4
%    semilogy(1:Nit, reshape(error(4,:), 1, Nit), 'kx:', 'linewidth', lw, 'markersize', ms);
%    legend('Gauss-Legendre','Gauss-Lobatto','Gauss-Radau','Equi-noleft','location','northeast'); legend boxoff;
% else
%     legend('Gauss-Legendre','Gauss-Lobatto','Gauss-Radau','location','northeast'); legend boxoff;
% end
% ylim([1e-16 1e2]);
% xlim([1 Nit]);
% ylabel('Iteration error','fontsize',fs);
% xlabel('Iteration','fontsize',fs);
% title(['N_x = ', num2str(Nx)]);
% set(gca,'fontsize',fs);
% 
fig = figure(2);
set(fig, 'Toolbar','none');
set(fig, 'PaperPositionMode', 'auto');
set(fig, 'Units', 'centimeters');
set(fig, 'OuterPosition', [0, 0, 16, 16]);
set(fig, 'Color','white');
semilogy(1:Nit, reshape(err_end(1,:), 1, Nit), 'ro:', 'linewidth', lw, 'markersize', ms); hold on;
semilogy(1:Nit, reshape(err_end(2,:), 1, Nit), 'bo:', 'linewidth', lw, 'markersize', ms);
legend('Gauss-Legendre','Gauss-Chebyshev', 'location','northeast'); legend boxoff;
%semilogy(1:Nit, reshape(err_end(3,:), 1, Nit), 'go:', 'linewidth', lw, 'markersize', ms);
% if nr_types==4
%    semilogy(1:Nit, reshape(err_end(4,:), 1, Nit), 'kx:', 'linewidth', lw, 'markersize', ms);
%    legend('Gauss-Legendre','Gauss-Lobatto','Gauss-Radau','Equi-noleft','location','southwest'); legend boxoff;
% else
%     legend('Gauss-Legendre','Gauss-Lobatto','Gauss-Radau','location','southwest'); legend boxoff;
% end
ylim([1e-16 1e2]);
xlim([1 Nit]);
ylabel('PDE error','fontsize',fs);
xlabel('Iteration','fontsize',fs);
title(['N_x = ', num2str(Nx)]);
set(gca,'fontsize',fs);
% 
% subplot(133)
% semilogy(1:Nit, reshape(residual(1,:), 1, Nit), 'ro:', 'linewidth', lw, 'markersize', ms); hold on;
% semilogy(1:Nit, reshape(residual(2,:), 1, Nit), 'bo:', 'linewidth', lw, 'markersize', ms);
% semilogy(1:Nit, reshape(residual(3,:), 1, Nit), 'go:', 'linewidth', lw, 'markersize', ms);
% if nr_types==4
%    semilogy(1:Nit, reshape(residual(4,:), 1, Nit), 'kx:', 'linewidth', lw, 'markersize', ms);
%    legend('Gauss-Legendre','Gauss-Lobatto','Gauss-Radau','Equi-noleft','location','northeast'); legend boxoff;
% else
%     legend('Gauss-Legendre','Gauss-Lobatto','Gauss-Radau','location','northeast'); legend boxoff;
% end
% ylim([1e-16 1e2]);
% xlim([1 Nit]);
% ylabel('Residual','fontsize',fs);
% xlabel('Iteration','fontsize',fs);
% title(['N_x = ', num2str(Nx)]);
% set(gca,'fontsize',fs);

% filename = ['sdc_lambda', num2str(lambda) '_M', num2str(nodes) '_.eps'];
% print('-depsc', filename);