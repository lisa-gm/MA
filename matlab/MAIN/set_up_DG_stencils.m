%% setting up the DG matrices for time stepping

function [K_t, M_t, N_t] = set_up_DG_stencils(dt)
%% time operators
syms t

%% basis function on reference element q = 1 
l1(t) = 1-t/dt;
l2(t) = t/dt;
basis = {l1,l2};

K_t = zeros(2);
M_t = zeros(2);
N_t = zeros(2);

for k = 1:2
    for l = 1:2
        
        bk = basis{k};
        bl = basis{l};
        
        K_t(k,l) = - int(diff(bk,t)*bl,0,dt) + bk(dt)*bl(dt);              
        
        M_t(k,l) = int(bk*bl,0,dt);      
        N_t(k,l) = bk(0)*bl(dt);       
    end
end