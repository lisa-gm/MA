clear all;
close all;
beep off;

addpath('../src');
addpath('../src/quad_nodes');
addpath('../../mparareal_oo/src');

Tend   = 1.0;
nodes  = 15;
lambda = -5.0 + 0i;

coll    = collocation_sdc(0, Tend, nodes, 'gauss-legendre');

yrhs = ones(nodes,1);

Nit = 50;

Q_pic = lambda*coll.Qmat;
Q_sdc_pre = ( speye(nodes) - lambda*coll.Qdelta_mat );
Q_sdc_left = ( lambda*(coll.Qmat - coll.Qdelta_mat) );

[ef_pic, e_pic] = eig(Q_pic);
[ef_sdc, e_sdc] = eig(Q_sdc_pre\Q_sdc_left);

[sc_pic ind_pic] = max(abs(diag(e_pic)));
[sc_sdc ind_sdc] = max(abs(diag(e_sdc)));

yex = (speye(nodes) - Q_pic)\yrhs;

y0   = ones(nodes,1);
%y0   = yex;
y_sdc = y0;
y_pic = y0;
err   = zeros(2,Nit);

coeff_0_pic = ef_pic\(y0 - yex);
coeff_0_sdc = ef_sdc\(y0 - yex);

coeff_sdc = zeros(nodes,Nit);
coeff_pic = zeros(nodes,Nit);

for nn=1:Nit
    y_sdc = Q_sdc_pre\( yrhs + Q_sdc_left*y_sdc);
    y_pic = Q_pic*y_pic + yrhs;
    
    err_sdc = y_sdc - yex;
    err_pic = y_pic - yex;
    

    err(1,nn) = norm(err_sdc, inf);
    err(2,nn) = norm(err_pic, inf);
    
    % Now project error to eigenspaces
    coeff_sdc(1:nodes,nn) = ef_sdc\err_sdc;
    coeff_pic(1:nodes,nn) = ef_pic\err_pic;
    
    
%     fig = figure(1); clf;
%     set(fig, 'Toolbar','none');
%     set(fig, 'PaperPositionMode', 'auto');
%     set(fig, 'Units', 'centimeters');
%     set(fig, 'OuterPosition', [0, 0, 16, 16]);
%     set(fig, 'Color','white');
%     semilogy(1:nodes, abs(coeff_sdc(:,nn)), 'ro--', 1:nodes, abs(coeff_pic(:,nn)), 'bd--');
%     ylim([1e-16 1e6]);
%     xlim([1 nodes]);
%     drawnow;
%     pause(0.8);
    
end

fig = figure(2);
set(fig, 'Toolbar','none');
set(fig, 'PaperPositionMode', 'auto');
set(fig, 'Units', 'centimeters');
set(fig, 'OuterPosition', [0, 0, 16, 16]);
set(fig, 'Color','white');

semilogy(1:Nit, abs(coeff_sdc)); hold on;
semilogy(1:Nit, max(abs(coeff_sdc(:,1)))*sc_sdc.^(0:Nit-1), 'k--', 'linewidth', 2);
semilogy(1:Nit, err(1,:), 'r--', 'linewidth', 2);
semilogy(1:Nit, sum(abs(coeff_sdc),1), 'g:','linewidth',2);
ylim([1e-18 1e15]);

fig = figure(3);
set(fig, 'Toolbar','none');
set(fig, 'PaperPositionMode', 'auto');
set(fig, 'Units', 'centimeters');
set(fig, 'OuterPosition', [0, 0, 16, 16]);
set(fig, 'Color','white');

semilogy(1:Nit, abs(coeff_pic)); hold on;
semilogy(1:Nit, max(abs(coeff_pic(:,1)))*sc_pic.^(0:Nit-1), 'k--', 'linewidth', 2);
semilogy(1:Nit, err(2,:), 'r--', 'linewidth', 2);
semilogy(1:Nit, sum(abs(coeff_pic),1), 'g:','linewidth',2);
ylim([1e-18 1e15]);

fig = figure(4);
set(fig, 'Toolbar','none');
set(fig, 'PaperPositionMode', 'auto');
set(fig, 'Units', 'centimeters');
set(fig, 'OuterPosition', [0, 0, 16, 16]);
set(fig, 'Color','white');

ev_pic_norm = zeros(1,nodes);
ev_sdc_norm = zeros(1,nodes);
for nn=1:nodes
    ev_pic_norm(1,nn) = norm( coeff_0_pic(nn)*ef_pic(:,nn), inf);
    ev_sdc_norm(1,nn) = norm( coeff_0_sdc(nn)*ef_sdc(:,nn), inf);
end


semilogy(1:nodes, ev_sdc_norm, 'ro--', 1:nodes, ev_pic_norm, 'bd--'); hold on;
semilogy(1:nodes, zeros(1,nodes)+cond(inv(ef_pic), inf), 'b', 1:nodes, zeros(1,nodes)+cond(inv(ef_sdc), inf), 'r');
legend('SDC','Picard','location','southeast');
% The dot marks the eigenvector belonging to the maximum absolute
% eigenvalue (=spectral radius)
semilogy(ind_pic, ev_pic_norm(ind_pic), 'bd', 'markerfacecolor', 'b', 'markersize', 10);
semilogy(ind_sdc, ev_sdc_norm(ind_sdc), 'ro', 'markerfacecolor', 'r', 'markersize', 10);
xlabel('m');
ylabel('abs(a_m)');
%ylim([1e-18 1e15]);
