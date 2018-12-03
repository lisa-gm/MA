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
addpath('../../matrix_mg');

Nx = 501;
N_half = (Nx+1)/2;

nu = 1.0;
[A, mesh] = getLapMat(0, 1, Nx);
A         = nu*A;

y0        = sin(pi*mesh).';
y_ex_fh   = @(t) exp(-nu*pi^2*t)*y0;

Tend   = 0.25;
Nsteps = 1;
Tmesh  = linspace(0,Tend,Nsteps+1);
nodes  = 5;
Nit    = 25;
method = {'jac', 'full'};
% method = {'jac-nodamp', 'jac-nodamp'};

sol_end = solution_arraylinear(y_ex_fh(Tend), A, speye(Nx));

sol_0_a   = solution_arraylinear_iter(y0,       A, speye(Nx), 0, 1, nu, method, 1);
sol_0_b   = solution_arraylinear_iter(y0,       A, speye(Nx), 0, 1, nu, {'full', 'full'}, 1);

colls{1}    = collocation_sdc(0, Tend, nodes, 'gauss-lobatto');
colls{1}    = colls{1}.initialize(sol_0_a);

colls{2} = collocation_sdc(0, Tend, nodes, 'gauss-lobatto');
colls{2} = colls{2}.initialize(sol_0_b);

nr_types = length(colls);

% For each set of nodes, generate cell array that contains the exact
% solution at the nodes
sol_ex      = cell(nr_types,nodes);
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

err      = zeros(Nit, nr_types);
err_spec = zeros(N_half, Nit, nr_types);

% Now perfom Nit iterations
for nn=1:Nit
    for ii=1:nr_types
        
        % Do Picard iteration
        colls{ii} = colls{ii}.applySweep(1.0);
        
        % Fetch approximation at Tend
        uend            = colls{ii}.get_end_value;
        
        % Compute difference to analytical solution and compute error at end
        diff       = uend.axpy(-1.0, sol_end);
        err(nn,ii) = diff.getnorm;
        
        diff_spec = abs(fftshift(fft(diff.y)));
        err_spec(:,nn,ii) = fliplr(diff_spec(1:N_half));
        
    end
    
end

fig = figure(1);
set(fig, 'Toolbar','none');
set(fig, 'PaperPositionMode', 'auto');
set(fig, 'Units', 'centimeters');
set(fig, 'OuterPosition', [0, 0, 16, 16]);
set(fig, 'Color','white');

offset = Nit-1;

for nn=1:(Nit-offset)
    clf;
    semilogy(1:N_half, err_spec(:,nn+offset,1), 'r', 'linewidth', lw); hold on;
    semilogy(1:N_half, err_spec(:,nn+offset,2), 'k', 'linewidth', lw);
    legend('ISDC','SDC','location','southeast');
    e = zeros(1,N_half);
    semilogy(1:N_half, e+err(nn+offset,1), 'r--');
    semilogy(1:N_half, e+err(nn+offset,2), 'k--');
    ylim([1e-10 1e5]);
    title(['k=',num2str(nn+offset)]);
    xlabel('Wave number', 'fontsize', fs);
    ylabel('abs(a_n)', 'fontsize', fs);
    set(gca,'fontsize',fs);
    drawnow;
    
    if nn<10
        filename = ['figs/isdc_errspec_0', num2str(nn) 'png'];
    else
        filename = ['figs/isdc_errspec_', num2str(nn) 'png'];
    end
    print('-dpng',filename);
    
end