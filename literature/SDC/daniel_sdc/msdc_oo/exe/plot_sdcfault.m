clear all;
close all;
beep off;

lw = 1.0;
fs = 12;
ms = 10;

addpath('../src');
addpath('../src/quad_nodes');
addpath('../../mparareal_oo/src');

% Equation (3.4), page 15
lambda  = 1.0;%*1i;
f       = @(y) lambda*y;
Tend    = 1.0;
y0      = 1.0;

% Perturbation
s = 1.5;

% Number of iterations
Nit = 20;

y_ex_fh = @(t) y0*exp(lambda*t);

nodes        = 3;
q_norm       = zeros(3,1);
q_delta_norm = zeros(3,1);

y_0     = y0*ones(nodes,1);
sol_0   = solution_arraylinear(y0, lambda, 1.0);

%colls    = collocation_sdc(0, Tend, nodes, 'gauss-legendre');
colls    = collocation_sdc(0, Tend, nodes, 'gauss-lobatto');
%colls    = collocation_sdc(0, Tend, nodes, 'gauss-radau');
colls    = colls.initialize(sol_0);

lambda_mat = lambda*speye(nodes);

% For the linear test equation y' = lambda*y, a sweep of SDC(theta,gamma)
% can be easily written as a matrix so that u_np1 = Mit*u_n0 + M0*u0
Mit = ( speye(nodes) - colls.Qdelta_mat*lambda_mat ) \ ( ( colls.Qmat - colls.Qdelta_mat)*lambda_mat );
M0  = ( speye(nodes) - colls.Qdelta_mat*lambda_mat ) \ ( speye(nodes) );

% Verify correctness of composed iteration matrices
%[Mit_ref, M0_ref] = colls.getScalarSweepMatrix(lambda, 1.0, 1.0);
% norm(Mit - Mit_ref, inf)
% norm(M0 - M0_ref, inf)

s_vec = logspace(-1, 4, 30);
err_final = zeros(3,length(s_vec));

for ii=1:length(s_vec)
    
    s = s_vec(ii);
    
    % Introduce perturbation of lambda at second node
    lambda_mat_pert = lambda*speye(nodes);
    lambda_mat_pert(2,2) = s*lambda;
    Mit_pert = ( speye(nodes) - colls.Qdelta_mat*lambda_mat_pert ) \ ( ( colls.Qmat - colls.Qdelta_mat)*lambda_mat_pert );
    M0_pert  = ( speye(nodes) - colls.Qdelta_mat*lambda_mat_pert ) \ ( speye(nodes) );
    
    yy = y_0;
    yy_pert = y_0;
    for iter=1:Nit
        yy = Mit*yy + M0*y_0;
        
        % apply perturbation in third iteration only
        if (iter==3)
            yy_pert = Mit_pert*yy_pert + M0_pert*y_0;
        else
            yy_pert = Mit*yy_pert + M0*y_0;
        end
        
        if (iter==3)
            err_final(1,ii) = abs(yy_pert(end) - y_ex_fh(Tend));
        elseif (iter==10)
            err_final(2,ii) = abs(yy_pert(end) - y_ex_fh(Tend));
                        
        elseif (iter==20)
            err_final(3,ii) = abs(yy_pert(end) - y_ex_fh(Tend));
        end
        
    end
    
    
    
end

fig = figure(1);
set(fig, 'Toolbar','none');
set(fig, 'PaperPositionMode', 'auto');
set(fig, 'Units', 'centimeters');
set(fig, 'OuterPosition', [0, 0, 32, 32]);
set(fig, 'Color','white');


specrad = max(abs(eig(full(Mit))));

loglog(s_vec, err_final(1,:), 'bo', 'markersize', ms, 'linewidth', lw, 'markerfacecolor', 'b'); hold on;
loglog(s_vec, err_final(2,:), 'r>', 'markersize', ms, 'linewidth', lw, 'markerfacecolor', 'r'); 
loglog(s_vec, err_final(3,:), 'gs', 'markersize', ms, 'linewidth', lw, 'markerfacecolor', 'g'); 
loglog(s_vec, err_final(1,:)*specrad^(10-3), 'r-', 'linewidth', lw); 
loglog(s_vec, err_final(1,:)*specrad^(20-3), 'g-', 'linewidth', lw);
legend('Sweep 3','Sweep 10', 'Sweep 20', 'Theory sweep 10', 'Theory sweep 20','location','NorthWest'); legend boxoff;
set(gca, 'fontsize', fs);
%
%
% tt = linspace(0, 1, 100);
% plot(colls.nodes, yy, 'bx-'); hold on;
% plot(colls.nodes, yy_pert, 'ro-', 'markerfacecolor' ,'r');
% plot(tt, y_ex_fh(tt), 'b', 'linewidth', 1.5);