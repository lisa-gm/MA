clear all;
close all;

addpath('../src');
addpath('../exe');

T0 = 0.0;
T1 = 0.5;
nodes = {'gauss-legendre', 'gauss-lobatto', 'gauss-radau'};

% First test: Verify that using identical set of nodes on both levels
% returns the identity
for tt=1:length(nodes)
    for M=2:8
        coll_f = collocation(T0, T1, M, nodes{tt});
        [R, P] = get_pr(coll_f.nodes, coll_f.nodes);
        assert( norm(R - speye(M), inf) < 1e-13, 'R is not the identity');
        assert( norm(P - speye(M), inf) < 1e-13, 'P is not the identity');
    end
end

% Define function for tests below
%u_fh   = @(t) erf( (t - 0.5*(T1-T0))./0.1 );
u_fh = @(t) sin(3*pi*t).^2;

% Restriction and interpolation error versus M
Mmax = 20;
err  = zeros(2, Mmax-1);
for m=2:Mmax
    coll_f = collocation(T0, T1, m+1, 'gauss-legendre');
    coll_c = collocation(T0, T1, m, 'gauss-radau');
    [R, P] = get_pr(coll_f.nodes, coll_c.nodes);
    
    % compute prolongated solution
    U_c   = u_fh(coll_c.nodes);
    U_c_p = P*U_c;
    
    % compute restricted solution
    U_f   = u_fh(coll_f.nodes);
    U_f_r = R*U_f;
    
    err(1,m-1) = norm(U_c - U_f_r, inf);
    err(2,m-1) = norm(U_f - U_c_p, inf);
end

figure(1);
semilogy(2:Mmax, err(1,:), 'bo-', 2:Mmax, err(2,:), 'rs-');

% Eye norm test
% coll_f = collocation(T0, T1, 12, 'gauss-legendre');
% coll_c = collocation(T0, T1, 6, 'gauss-radau');
% [R, P] = get_pr(coll_f.nodes, coll_c.nodes);
% 
% U_c   = u_fh(coll_c.nodes);
% U_c_p = P*U_c;
% 
% U_f   = u_fh(coll_f.nodes);
% U_f_r = R*U_f;
% 
% figure(2)
% subplot(211)
% plot(coll_f.nodes, U_f, 'bo', coll_f.nodes, U_c_p, 'r--');
% legend('Fine solution','Interpolated'); legend boxoff;
% 
% subplot(212)
% plot(coll_c.nodes, U_c, 'bo', coll_c.nodes, U_f_r, 'r--');
% legend('Coarse solution', 'Restricted'); legend boxoff;
% 
% fprintf('Fine level error: %5.3e \n', norm(U_f - U_c_p, inf));
% fprintf('Coarse level error: %5.3e \n', norm(U_c - U_f_r, inf));