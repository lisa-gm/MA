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
re_v = linspace(-10, 0, 110);
im_v = linspace(0, 5,    60);

Tend    = 1.0;

nodes    = 3;

coll    = collocation_sdc(0, Tend, nodes, 'gauss-legendre');
bound   = zeros(length(im_v), length(re_v));

for jj=1:length(im_v)
    
    for ii=1:length(re_v)
        
        lambda = re_v(ii) + 1i*im_v(jj);
        
%         temp = 0;
%         for mm=1:nodes
%               temp = temp + 2./abs(1 - lambda*coll.delta_m(mm));
%         end
        M = (speye(nodes) - lambda*coll.Qdelta_mat)\speye(nodes);
        
        bound(jj,ii) = norm(M, inf);
        
    end
    
end

fig = figure(1); clf;
set(fig, 'Toolbar','none');
set(fig, 'PaperPositionMode', 'auto');
set(fig, 'Units', 'centimeters');
set(fig, 'OuterPosition', [0, 0, 16, 16]);
set(fig, 'Color','white');

surf(re_v, im_v, bound); colorbar;
view(0,90);