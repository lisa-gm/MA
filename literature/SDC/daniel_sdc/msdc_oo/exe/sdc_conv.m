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
addpath('~/Programs/Matlab/ExampleCodes/LinearParareal');

% Solves y' = lambda*y with Picard iteration
lambda  = -1.0 + 0.0*1i;

mass    = speye(size(lambda)); % probably can't be a different value
f       = @(y) lambda*y;

Tend    = 1.0;
y0      = 1.0;
Nit_max = 5;
nodes   = 5;

y_ex_fh = @(t) y0*expm(lambda*t);

sol_0   = solution_arraylinear(y0,            lambda, mass);
sol_ex  = solution_arraylinear(y_ex_fh(Tend), lambda, mass);

Nsteps_v = 5:10;

nr_types = 4;
err_end  = zeros(nr_types, length(Nsteps_v));

for ll=1:length(Nsteps_v)
    
    Nsteps = Nsteps_v(ll);
    Tmesh  = linspace(0,Tend,Nsteps+1);
    
    colls       = cell(nr_types,1);
    colls{1}    = summed_collocation(Tmesh, nodes, 'gauss-legendre');
    colls{2}    = summed_collocation(Tmesh, nodes, 'gauss-lobatto');
    colls{3}    = summed_collocation(Tmesh, nodes, 'gauss-radau');
    colls{4}    = summed_collocation(Tmesh, nodes, 'equi-noleft');
    
    theta       = 1.0*ones(1,nr_types);
          
    gamma       = 1.0;
    
    % Generate cell array with sol_0 at all nodes as initial value
    y_0   = cell(1,nodes);
    for nn=1:nodes
        y_0{1,nn} = sol_0;
    end
    
    for ii=1:nr_types
        colls{ii} = colls{ii}.initialize(sol_0,1);
    end
    
    
    % Now perfom Nit iterations
    for ii=1:nr_types
        for nn=1:Nsteps
            
            res = 1.0;
            nit = 0;
            while nit<=Nit_max
                % Do sweeps
                colls{ii} = colls{ii}.applySweep(theta(ii), gamma, nn);
                res = colls{ii}.collocations{nn}.get_residual;
                nit = nit+1;
            end
            fprintf('%s -- #Iterations: %2i -- Residual: %5.2e \n', colls{ii}.collocations{1}.type, nit, res);
            
            if nn<Nsteps
                
                % Fetch final solution from current timestep
                sol_end = colls{ii}.collocations{nn}.get_end_value;
                
                % and initialize next timestep with new starting value; also updates node solutions
                colls{ii} = colls{ii}.initialize(sol_end, nn+1);
               
            end
        end
        
        % Get solution at TEnd
        sol_end = colls{ii}.collocations{Nsteps}.get_end_value;
        
        % difference to exact solution
        diff = sol_end.axpy(-1.0, sol_ex);
        
        % error
        err_end(ii,ll) = diff.getnorm;%/sol_ex.getnorm;
    end
    
end

fig = figure(1);
set(fig, 'Toolbar','none');
set(fig, 'PaperPositionMode', 'auto');
set(fig, 'Units', 'centimeters');
set(fig, 'OuterPosition', [0, 0, 16, 16]);
set(fig, 'Color','white');

loglog(Nsteps_v, err_end(1,:), 'ro:', 'linewidth', lw, 'markersize', ms); hold on;
loglog(Nsteps_v, err_end(2,:), 'bo:', 'linewidth', lw, 'markersize', ms);
loglog(Nsteps_v, err_end(3,:), 'go:', 'linewidth', lw, 'markersize', ms);
loglog(Nsteps_v, err_end(4,:), 'kx:', 'linewidth', lw, 'markersize', ms);
legend('Gauss-Legendre','Gauss-Lobatto','Gauss-Radau','Equi-noleft','location','northeast'); legend boxoff;

order = (Nsteps_v(1)./Nsteps_v).^(Nit_max+1);
max_err_0 = max( max(err_end(1,1), err_end(2,1)), max(err_end(3,1),err_end(4,1)) );
loglog(Nsteps_v, 1.25*max_err_0*order, 'k-', 'linewidth', lw);

xlim([min(Nsteps_v) max(Nsteps_v)]);
xlabel('timesteps','fontsize', fs);
ylabel('Relative error', 'fontsize', fs);
title(['\lambda=', num2str(lambda) ', M=',num2str(nodes)]);
set(gca,'fontsize',fs);