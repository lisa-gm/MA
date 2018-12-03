clear all;
close all;
beep off;

addpath('../src');

M = 15;
% type = 'gauss-lobatto';
type = 'gauss-legendre';
%type = 'equi';
% type = 'gauss-chebyshev';

coll = collocation(-1, 1, M, type);
nodes = coll.nodes;

xplot = linspace(-1,1,200);
yplot = zeros(M, length(xplot));

for mm=1:M
    ynodes = zeros(1,length(nodes));
    ynodes(1,mm) = 1.0;
    
    lagr = interp_poly(nodes, ynodes);
    yplot(mm,:) = lagr.evaluate(xplot);
end

plot(xplot, abs(yplot)); hold on;
plot(xplot, sum(yplot, 1), 'k-', 'linewidth', 1.5);