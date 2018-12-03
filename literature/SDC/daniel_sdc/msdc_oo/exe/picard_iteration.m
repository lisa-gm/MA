clear all;
close all;
beep off;

lw = 1.4;
fs = 12;
ms = 12;

addpath('../src');
addpath('../src/quad_nodes');
addpath('../../mparareal_oo/src');

% Solves y' = lambda*y with Picard iteration
lambda  = -1.0 + 0.0*1i;
mass    = 1.0; % probably can't be a different value
f       = @(y) lambda*y;
Tend    = 1.0;
y0      = 1.0;

y_ex_fh = @(t) y0*exp(lambda*t);

theta  = [0.0 1.0];

Nsteps = 1;
Tmesh  = linspace(0,Tend,Nsteps+1);
nodes  = 9;
Nit    = 25;

sol_0   = solution_arraylinear(y0, lambda, mass);
sol_end = solution_arraylinear(y_ex_fh(Tend), lambda, mass);

%tt = theta values (hardcoded to be two here)
eigv = zeros(2,3);
for tt=1:2
    colls{tt,1}    = collocation_sdc(0, Tend, nodes, 'gauss-legendre');
    colls{tt,1}    = colls{tt,1}.initialize(sol_0);
    Mtemp          = ( speye(nodes) - lambda*theta(tt)*colls{tt,1}.Qdelta_mat )\( lambda*(colls{tt,1}.Qmat - theta(tt)*colls{tt,1}.Qdelta_mat) );
    eigv(tt,1)     = max(abs(eig(full(Mtemp))));
    %eigv(tt,1)     = svds(Mtemp,1);
    
    colls{tt,2}    = collocation_sdc(0, Tend, nodes, 'gauss-lobatto');
    colls{tt,2}    = colls{tt,2}.initialize(sol_0);
    Mtemp          = ( speye(nodes) - lambda*theta(tt)*colls{tt,2}.Qdelta_mat )\( lambda*(colls{tt,2}.Qmat - theta(tt)*colls{tt,2}.Qdelta_mat) );
    eigv(tt,2)     = max(abs(eig(full(Mtemp))));
    %eigv(tt,2)     = svds(Mtemp,1);
    
    colls{tt,3}    = collocation_sdc(0, Tend, nodes, 'gauss-radau');
    colls{tt,3}    = colls{tt,3}.initialize(sol_0);
    Mtemp          = ( speye(nodes) - lambda*theta(tt)*colls{tt,3}.Qdelta_mat )\( lambda*(colls{tt,3}.Qmat - theta(tt)*colls{tt,3}.Qdelta_mat) );
    eigv(tt,3)     = max(abs(eig(full(Mtemp))));
    %eigv(tt,3)     = svds(Mtemp,1);

end

% For each set of nodes, generate cell array that contains the exact
% solution at the nodes
sol_ex      = cell(3,nodes);
for ii=1:3
    for nn=1:nodes
        sol_ex{ii,nn} = solution_arraylinear(y_ex_fh(colls{1,ii}.nodes(nn)), lambda, mass);
    end
end

% For each set of nodes, compute also the correct collocation solution by
% solving U - Q*lambda*U = U0
sol_coll_ex = cell(3,nodes);
for ii=1:3
    u_coll = (speye(nodes) - colls{1,ii}.Qmat*lambda)\(y0*ones(nodes,1));
    for nn=1:nodes
       sol_coll_ex{ii,nn} = solution_arraylinear( u_coll(nn,1), lambda, mass);
    end
end

% Generate cell array with sol_0 at all nodes as initial value
y_0   = cell(1,nodes);
for nn=1:nodes
    y_0{1,nn} = sol_0;
    %y_0{1,nn} = sol_coll_ex{1,nn};
end
for tt=1:2
    for ii=1:3
        colls{tt,ii} = colls{tt,ii}.set_node_solutions(y_0);
    end
end

error    = zeros(2,3,Nit);
err_end  = zeros(2,3,Nit);
residual = zeros(2,3,Nit);

% Now perfom Nit iterations
for nn=1:Nit
    for tt=1:2
        for ii=1:3
            
            % Do Picard iteration
            colls{tt,ii} = colls{tt,ii}.applySweep(theta(tt));
            
            % Fetch approximation at Tend
            uend            = colls{tt,ii}.get_end_value;
            
            % Compute difference to analytical solution and compute error
            diff            = uend.axpy(-1.0, sol_end);
            err_end(tt,ii,nn)  = diff.getnorm;
            
            % Fetch residual
            residual(tt,ii,nn) = colls{tt,ii}.get_residual;
            
            % compute defect on all nodes to analytical solution
            def = zeros(nodes,1);
            for mm=1:nodes
%                 diff      = sol_ex{ii,mm}.axpy(-1.0,
%                 colls{tt,ii}.node_solutions{1,mm});
                diff      = sol_coll_ex{ii,mm}.axpy(-1.0, colls{tt,ii}.node_solutions{1,mm});
                def(mm,1) = diff.getnorm;
            end
            error(tt,ii,nn) = norm(def, inf);
        end
    end
end

fig = figure(1);
set(fig, 'Toolbar','none');
set(fig, 'PaperPositionMode', 'auto');
set(fig, 'Units', 'centimeters');
set(fig, 'OuterPosition', [0, 0, 36, 24]);
set(fig, 'Color','white');

%subplot(231)
subplot(121)
semilogy(1:Nit, reshape(error(1,1,:), 1, Nit), 'ro:', 'linewidth', lw, 'markersize', ms); hold on;
semilogy(1:Nit, reshape(error(1,2,:), 1, Nit), 'bo:', 'linewidth', lw, 'markersize', ms);
semilogy(1:Nit, reshape(error(1,3,:), 1, Nit), 'go:', 'linewidth', lw, 'markersize', ms);
legend('Gauss-Legendre','Gauss-Lobatto','Gauss-Radau','location','northeast'); legend boxoff;
semilogy(1:Nit, error(1,1,end)*eigv(1,1).^(0:Nit-1)/eigv(1,3)^(Nit-1), 'r-', 'linewidth', lw);
semilogy(1:Nit, error(1,2,end)*eigv(1,2).^(0:Nit-1)/eigv(1,2)^(Nit-1), 'b-', 'linewidth', lw);
semilogy(1:Nit, error(1,3,end)*eigv(1,3).^(0:Nit-1)/eigv(1,3)^(Nit-1), 'g-', 'linewidth', lw);
ylim([1e-16 1e0]);
xlim([1 Nit]);
ylabel('Iteration error','fontsize',fs);
xlabel('Iteration','fontsize',fs);
title(['\Theta = ', num2str(theta(1))], 'fontsize', fs);

% subplot(232)
% semilogy(1:Nit, reshape(err_end(1,1,:), 1, Nit), 'ro-', 'linewidth', lw, 'markersize', ms); hold on;
% semilogy(1:Nit, reshape(err_end(1,2,:), 1, Nit), 'bo-', 'linewidth', lw, 'markersize', ms);
% semilogy(1:Nit, reshape(err_end(1,3,:), 1, Nit), 'go-', 'linewidth', lw, 'markersize', ms);
% legend('Gauss-Legendre','Gauss-Lobatto','Gauss-Radau','location','northeast'); legend boxoff;
% ylabel('End value error','fontsize',fs);
% xlabel('Iteration','fontsize',fs);
% title(['\Theta = ', num2str(theta(1))], 'fontsize', fs);
% 
% subplot(233)
% semilogy(1:Nit, reshape(residual(1,1,:), 1, Nit), 'ro-', 'linewidth', lw, 'markersize', ms); hold on;
% semilogy(1:Nit, reshape(residual(1,2,:), 1, Nit), 'bo-', 'linewidth', lw, 'markersize', ms);
% semilogy(1:Nit, reshape(residual(1,3,:), 1, Nit), 'go-', 'linewidth', lw, 'markersize', ms);
% legend('Gauss-Legendre','Gauss-Lobatto','Gauss-Radau','location','northeast'); legend boxoff;
% ylabel('Residual','fontsize',fs);
% xlabel('Iteration','fontsize',fs);
% title(['\Theta = ', num2str(theta(1))], 'fontsize', fs);
%

% subplot(234)
subplot(122)
semilogy(1:Nit, reshape(error(2,1,:), 1, Nit), 'ro:', 'linewidth', lw, 'markersize', ms); hold on;
semilogy(1:Nit, reshape(error(2,2,:), 1, Nit), 'bo:', 'linewidth', lw, 'markersize', ms);
semilogy(1:Nit, reshape(error(2,3,:), 1, Nit), 'go:', 'linewidth', lw, 'markersize', ms);
legend('Gauss-Legendre','Gauss-Lobatto','Gauss-Radau','location','northeast'); legend boxoff;
semilogy(1:Nit, error(2,1,end)*eigv(2,1).^(0:Nit-1)/eigv(2,1)^(Nit-1), 'r-', 'linewidth', lw);
semilogy(1:Nit, error(2,2,end)*eigv(2,2).^(0:Nit-1)/eigv(2,2)^(Nit-1), 'b-', 'linewidth', lw);
semilogy(1:Nit, error(2,3,end)*eigv(2,3).^(0:Nit-1)/eigv(2,3)^(Nit-1), 'g-', 'linewidth', lw);
ylim([1e-16 1e0]);
xlim([1 Nit]);
ylabel('Iteration error','fontsize',fs);
xlabel('Iteration','fontsize',fs);
title(['\Theta = ', num2str(theta(2))], 'fontsize', fs);

% subplot(235)
% semilogy(1:Nit, reshape(err_end(2,1,:), 1, Nit), 'ro-', 'linewidth', lw, 'markersize', ms); hold on;
% semilogy(1:Nit, reshape(err_end(2,2,:), 1, Nit), 'bo-', 'linewidth', lw, 'markersize', ms);
% semilogy(1:Nit, reshape(err_end(2,3,:), 1, Nit), 'go-', 'linewidth', lw, 'markersize', ms);
% legend('Gauss-Legendre','Gauss-Lobatto','Gauss-Radau','location','northeast'); legend boxoff;
% ylabel('End value error','fontsize',fs);
% xlabel('Iteration','fontsize',fs);
% title(['\Theta = ', num2str(theta(2))], 'fontsize', fs);
% 
% subplot(236)
% semilogy(1:Nit, reshape(residual(2,1,:), 1, Nit), 'ro-', 'linewidth', lw, 'markersize', ms); hold on;
% semilogy(1:Nit, reshape(residual(2,2,:), 1, Nit), 'bo-', 'linewidth', lw, 'markersize', ms);
% semilogy(1:Nit, reshape(residual(2,3,:), 1, Nit), 'go-', 'linewidth', lw, 'markersize', ms);
% legend('Gauss-Legendre','Gauss-Lobatto','Gauss-Radau','location','northeast'); legend boxoff;
% ylabel('Residual','fontsize',fs);
% xlabel('Iteration','fontsize',fs);
% title(['\Theta = ', num2str(theta(2))], 'fontsize', fs);