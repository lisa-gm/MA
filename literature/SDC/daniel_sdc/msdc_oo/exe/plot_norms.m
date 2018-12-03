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
lambda  = -10 + 0*1i;
f       = @(y) lambda*y;
Tend    = 1.0;
y0      = 1.0;

y_ex_fh = @(t) y0*exp(lambda*t);

nodes_v      = 25;
q_norm       = zeros(3, length(nodes_v));
q_delta_norm = zeros(3, length(nodes_v));

estimate     = zeros(3, length(nodes_v),2);

theta = 0.5;
gamma = 1.0;
        
for m=1:length(nodes_v)
    nodes = nodes_v(m);
    
    Nsteps = 1;
    Tmesh  = linspace(0,Tend,Nsteps+1);
    Nit    = 20;
    
    y_0     = y0*ones(nodes,1);
    sol_0   = solution_arraylinear(y0, lambda, 1.0);
    
    colls{1}    = collocation_sdc(0, Tend, nodes, 'gauss-legendre');
    colls{1}    = colls{1}.initialize(sol_0);
        
    T = tril(ones(nodes_v(m)));
    
    for ii=1:1
        
        Q_pic = lambda*colls{ii}.Qmat;
        Q_sdc = colls{1}.getScalarSweepMatrix(lambda, theta, gamma);

        Q_other = ( speye(nodes) - theta*lambda*colls{ii}.Qdelta_mat )\( (1 - gamma)*speye(nodes) + lambda*(colls{ii}.Qmat - theta*colls{ii}.Qdelta_mat) );
                
        fig = figure(1);
        set(fig, 'Toolbar','none');
        set(fig, 'PaperPositionMode', 'auto');
        set(fig, 'Units', 'centimeters');
        set(fig, 'OuterPosition', [0, 0, 16, 16]);
        set(fig, 'Color','white');
        
        [ef_pic, e_pic] = eig(Q_pic);
        [ef_sdc, e_sdc] = eig(Q_sdc);
        [ef_other, e_other] = eig(Q_other);
        e_pic = diag(e_pic);
        e_sdc = diag(e_sdc);
        e_other = diag(e_other);
        
        fprintf('Spectral radius Picard: %5.3f \n', max(abs(e_pic)));
        fprintf('Spectral radius SDC:    %5.3f \n', max(abs(e_sdc)));
        fprintf('Spectral radius other:  %5.3f \n', max(abs(e_other)));
        
        plot(real(e_pic), imag(e_pic), 'bo', 'markersize', 8, 'markerfacecolor','b'); hold on;
        plot(real(e_sdc), imag(e_sdc), 'rd', 'markersize', 8, 'markerfacecolor','r');
        plot(real(e_other), imag(e_other), 'gs', 'markersize', 8, 'markerfacecolor','g');
        legend('Picard','SDC','SDC(\theta)');
        ra_pic = max(abs(e_pic));
        r_pic  = linspace(-ra_pic,ra_pic,100);
        plot(r_pic, sqrt(ra_pic^2 - r_pic.^2), 'b--');
        plot(r_pic, -sqrt(ra_pic^2 - r_pic.^2), 'b--');
        
        ra_sdc = max(abs(e_sdc));
        r_sdc  = linspace(-ra_sdc, ra_sdc, 100);
        plot(r_sdc, sqrt(ra_sdc^2 - r_sdc.^2), 'r--');
        plot(r_sdc, -sqrt(ra_sdc^2 - r_sdc.^2), 'r--');
        
        ra_other = max(abs(e_other));
        r_other  = linspace(-ra_other, ra_other, 100);
        plot(r_other, sqrt(ra_other^2 - r_other.^2), 'g--');
        plot(r_other, -sqrt(ra_other^2 - r_other.^2), 'g--');
        
        %         r = linspace(-1,1,1000);
        %         plot(r, sqrt(1 - r.^2), 'k');
        %         plot(r, -sqrt(1 - r.^2), 'k');
        title(['abs(\lambda) * \delta t = ', num2str(abs(lambda*Tend)) ', M = ', num2str(nodes)]);
        
        %         mm_re_max = max( max(real(e_pic)), max(real(e_sdc)));
        %         mm_re_min = min( min(real(e_pic)), min(real(e_sdc)));
        %
        %         mm_im_max = max( max(imag(e_pic)), max(imag(e_sdc)));
        %         mm_im_min = min( min(imag(e_pic)), min(imag(e_sdc)));
        %
        %         xlim([ mm_re_min mm_re_max]);
        %         ylim([ mm_im_min mm_im_max]);
        
        max(abs(eig(Q_pic)))
        max(abs(eig(Q_sdc)))
        fprintf('*** \n');
        
        %         q_delta_norm(ii,m) = s_sdc;
        %
        %         assert( abs( norm(colls{ii}.Qmat - colls{ii}.Qdelta_mat) - s_sdc ) < 1e-14, 'SDC: Error in norm and singular value');
        %         assert( abs( s_sdc - norm( (colls{ii}.Qmat-colls{ii}.Qdelta_mat)*v_sdc) ) < 1e-14, 'SDC: Error in singular vector');
        
        q_norm(ii,m)       = norm(lambda*colls{ii}.Qmat, inf);
        %q_delta_norm(ii,m) = norm( (speye(nodes) - lambda*colls{ii}.Qdelta_mat)\(lambda*(colls{ii}.Qmat - colls{ii}.Qdelta_mat)), inf);
        q_delta_norm(ii,m) = norm( inv(speye(nodes) - lambda*colls{ii}.Qdelta_mat), inf);
        %    q_norm(ii,m)       = max(abs(eig(Q_pic)));
        %    q_delta_norm(ii,m) = max(abs(eig(Q_sdc)));
        
        %   normality_pic(ii,m) = norm( Q_pic.'*Q_pic - Q_pic*Q_pic.' )/norm(Q_pic.'*Q_pic);
        %   normality_sdc(ii,m) = norm( Q_sdc.'*Q_sdc - Q_sdc*Q_sdc.' )/norm(Q_sdc.'*Q_sdc);
        
        nodes_leftpoint = [0.0 ; colls{ii}.nodes.'];
        dtau            = diff(nodes_leftpoint);
        sum_dtau = sum(dtau);
        min_dtau = min(dtau);
        max_dtau = max(dtau);
        
        estimate(ii,m,1) = Tend*abs(lambda);
        % estimate(ii,m,1) = Tend*(2^nodes/(exp(1)*nodes*log10(nodes)));
        %estimate(ii,m,2) = Tend/nodes_v(m);
        estimate(ii,m,2) = (1 - abs(lambda)*Tend)^(-1)*abs(lambda)*Tend*2;
    end
end

% fig = figure(1);
% set(fig, 'Toolbar','none');
% set(fig, 'PaperPositionMode', 'auto');
% set(fig, 'Units', 'centimeters');
% set(fig, 'OuterPosition', [0, 0, 16, 16]);
% set(fig, 'Color','white');
% plot(colls{1}.nodes, v_pic, 'bo-'); hold on;
% plot(colls{1}.nodes, v_sdc, 'ro-');
%
% fig = figure(2);
% set(fig, 'Toolbar','none');
% set(fig, 'PaperPositionMode', 'auto');
% set(fig, 'Units', 'centimeters');
% set(fig, 'OuterPosition', [0, 0, 16, 16]);
% set(fig, 'Color','white');
% plot(colls{1}.nodes, u_pic, 'b'); hold on;
% plot(colls{1}.nodes, u_sdc, 'r');
%
% fig = figure(1);
% set(fig, 'Toolbar','none');
% set(fig, 'PaperPositionMode', 'auto');
% set(fig, 'Units', 'centimeters');
% set(fig, 'OuterPosition', [0, 0, 16, 16]);
% set(fig, 'Color','white');
% plot(nodes_v, q_norm(1,:), 'bo-', 'markerfacecolor', 'b'); hold on;
% %plot(nodes_v, q_norm(2,:), 'ro-', 'markerfacecolor', 'r'); hold on;
% %plot(nodes_v, q_norm(3,:), 'go-', 'markerfacecolor', 'g'); hold on;
% %legend('Legendre','Lobatto','Radau');
% plot(nodes_v, estimate(1,:,1), 'k--', 'linewidth', 2.0);
% xlabel('M','fontsize', fs);
% ylabel('||Q||', 'fontsize',fs);
% set(gca,'fontsize', fs);
% %
% fig = figure(2);
% set(fig, 'Toolbar','none');
% set(fig, 'PaperPositionMode', 'auto');
% set(fig, 'Units', 'centimeters');
% set(fig, 'OuterPosition', [0, 0, 16, 16]);
% set(fig, 'Color','white');
% plot(nodes_v, q_delta_norm(1,:), 'bo-', 'markerfacecolor', 'b'); hold on;
% %plot(nodes_v, q_delta_norm(2,:), 'ro-', 'markerfacecolor', 'r'); hold on;
% %plot(nodes_v, q_delta_norm(3,:), 'go-', 'markerfacecolor', 'g'); hold on;
% %legend('Legendre','Lobatto','Radau');
% plot(nodes_v,estimate(1,:,2), 'k--','linewidth',2.0);
% xlabel('M','fontsize', fs);
% ylabel('||Q - Q_{\Delta}||', 'fontsize',fs);
% set(gca,'fontsize', fs);


% fig = figure(1);
% set(fig, 'Toolbar','none');
% set(fig, 'PaperPositionMode', 'auto');
% set(fig, 'Units', 'centimeters');
% set(fig, 'OuterPosition', [0, 0, 16, 16]);
% set(fig, 'Color','white');

% plot(nodes_v, q_delta_norm(1,:), 'ro-', 'markerfacecolor', 'r'); hold on;
%plot(nodes_v, q_norm(1,:), 'bo-', 'markerfacecolor', 'b'); hold on; plot(nodes_v, q_norm(1,end)*(nodes_v(end)./(nodes_v)), 'k-');

%plot(nodes_v, normality_sdc(1,:), 'ro-', 'markerfacecolor', 'r'); hold on;
%plot(nodes_v, normality_pic(1,:), 'bo-', 'markerfacecolor', 'b');