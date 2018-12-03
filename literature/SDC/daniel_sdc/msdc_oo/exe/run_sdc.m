clear all;
close all;
beep off;

lw = 1.4;
fs = 12;
ms = 10;

addpath('../src');
addpath('../src/quad_nodes');
addpath('../../mparareal_oo/src');

% Solves y' = lambda*y with Picard iteration
lambda  = -1.0 + 0.0*1i;
mass    = 1.0; % probably can't be a different value
f       = @(y) lambda*y;
Tend    = 1.0;
y0      = 1.0;

y_ex_fh = @(t) y0*expm(lambda*t);

Nsteps = 1;
Tmesh  = linspace(0,Tend,Nsteps+1);
nodes  = 6;
Nit    = 25;

sol_0   = solution_arraylinear(y0,            lambda, mass);
sol_end = solution_arraylinear(y_ex_fh(Tend), lambda, mass);

eigv = zeros(1,3);

colls{1}    = collocation_sdc(0, Tend, nodes, 'gauss-legendre');
colls{2}    = collocation_sdc(0, Tend, nodes, 'gauss-lobatto');
colls{3}    = collocation_sdc(0, Tend, nodes, 'gauss-radau');
colls{4}    = collocation_sdc(0, Tend, nodes, 'equi-noleft');
nr_types = length(colls);

for ii=1:nr_types
   colls{ii} = colls{ii}.initialize(sol_0);
   Mtemp = ( speye(nodes) - lambda*colls{ii}.Qdelta_mat )\( lambda*(colls{ii}.Qmat - colls{ii}.Qdelta_mat) );
   eigv(ii) = max(abs(eig(full(Mtemp))));
end


% For each set of nodes, generate cell array that contains the exact
% solution at the nodes
sol_ex      = cell(3,nodes);
for ii=1:nr_types
    for nn=1:nodes
        sol_ex{ii,nn} = solution_arraylinear(y_ex_fh(colls{1,ii}.nodes(nn)), lambda, mass);
    end
end

% For each set of nodes, compute also the correct collocation solution by
% solving U - Q*lambda*U = U0
sol_coll_ex = cell(nr_types,nodes);
for ii=1:nr_types
    u_coll = (speye(nodes) - colls{1,ii}.Qmat*lambda)\(y0*ones(nodes,1));
    for nn=1:nodes
       sol_coll_ex{ii,nn} = solution_arraylinear( u_coll(nn,1), lambda, mass);
    end
end

% Generate cell array with sol_0 at all nodes as initial value
y_0   = cell(1,nodes);
for nn=1:nodes
    y_0{1,nn} = sol_0;
end

for ii=1:nr_types
    colls{ii} = colls{ii}.set_node_solutions(y_0);
end

error    = zeros(nr_types,Nit);
err_end  = zeros(nr_types,Nit);
residual = zeros(nr_types,Nit);

% Now perfom Nit iterations
for nn=1:Nit
    for ii=1:nr_types
        
        % Do Picard iteration
        colls{ii} = colls{ii}.applySweep(1.0, 1.0);
        
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

fig = figure(1);
set(fig, 'Toolbar','none');
set(fig, 'PaperPositionMode', 'auto');
set(fig, 'Units', 'centimeters');
set(fig, 'OuterPosition', [0, 0, 32, 24]);
set(fig, 'Color','white');

subplot(131)
semilogy(1:Nit, reshape(error(1,:), 1, Nit), 'ro:', 'linewidth', lw, 'markersize', ms); hold on;
semilogy(1:Nit, reshape(error(2,:), 1, Nit), 'bo:', 'linewidth', lw, 'markersize', ms);
semilogy(1:Nit, reshape(error(3,:), 1, Nit), 'go:', 'linewidth', lw, 'markersize', ms);
if nr_types==4
   semilogy(1:Nit, reshape(error(4,:), 1, Nit), 'kx:', 'linewidth', lw, 'markersize', ms);
   legend('Gauss-Legendre','Gauss-Lobatto','Gauss-Radau','Equi-noleft','location','northeast'); legend boxoff;
else
    legend('Gauss-Legendre','Gauss-Lobatto','Gauss-Radau','location','northeast'); legend boxoff;
end
ylim([1e-16 1e2]);
xlim([1 Nit]);
ylabel('Iteration error','fontsize',fs);
xlabel('Iteration','fontsize',fs);
title(['\lambda = ', num2str(lambda)]);
set(gca,'fontsize',fs);

subplot(132)
semilogy(1:Nit, reshape(err_end(1,:), 1, Nit), 'ro:', 'linewidth', lw, 'markersize', ms); hold on;
semilogy(1:Nit, reshape(err_end(2,:), 1, Nit), 'bo:', 'linewidth', lw, 'markersize', ms);
semilogy(1:Nit, reshape(err_end(3,:), 1, Nit), 'go:', 'linewidth', lw, 'markersize', ms);
if nr_types==4
   semilogy(1:Nit, reshape(err_end(4,:), 1, Nit), 'kx:', 'linewidth', lw, 'markersize', ms);
   legend('Gauss-Legendre','Gauss-Lobatto','Gauss-Radau','Equi-noleft','location','southwest'); legend boxoff;
else
    legend('Gauss-Legendre','Gauss-Lobatto','Gauss-Radau','location','southwest'); legend boxoff;
end
ylim([1e-16 1e2]);
xlim([1 Nit]);
ylabel('ODE error','fontsize',fs);
xlabel('Iteration','fontsize',fs);
title(['\lambda = ', num2str(lambda)]);
set(gca,'fontsize',fs);

subplot(133)
semilogy(1:Nit, reshape(residual(1,:), 1, Nit), 'ro:', 'linewidth', lw, 'markersize', ms); hold on;
semilogy(1:Nit, reshape(residual(2,:), 1, Nit), 'bo:', 'linewidth', lw, 'markersize', ms);
semilogy(1:Nit, reshape(residual(3,:), 1, Nit), 'go:', 'linewidth', lw, 'markersize', ms);
if nr_types==4
   semilogy(1:Nit, reshape(residual(4,:), 1, Nit), 'kx:', 'linewidth', lw, 'markersize', ms);
   legend('Gauss-Legendre','Gauss-Lobatto','Gauss-Radau','Equi-noleft','location','northeast'); legend boxoff;
else
    legend('Gauss-Legendre','Gauss-Lobatto','Gauss-Radau','location','northeast'); legend boxoff;
end
ylim([1e-16 1e2]);
xlim([1 Nit]);
ylabel('Residual','fontsize',fs);
xlabel('Iteration','fontsize',fs);
title(['\lambda = ', num2str(lambda)]);
set(gca,'fontsize',fs);

filename = ['sdc_lambda', num2str(lambda) '_M', num2str(nodes) '_.eps'];
print('-depsc', filename);