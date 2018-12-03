function [G, K, M] = set_up_FEM_1D_mat(N)
h = 1/(N-1);

%% time operators
syms t

%% basis function on reference element q = 1 
l1(t) = 1-t/h;
l2(t) = t/h;
basis = {l1,l2};

G_stencil = zeros(2);
K_stencil = zeros(2);
M_stencil = zeros(2);

for k = 1:2
    for l = 1:2
        
        bk = basis{k};
        bl = basis{l};
        
        G_stencil(k,l) = - int(diff(bk,t)*bl, 0,h);
        K_stencil(k,l) = int(diff(bk,t)*diff(bl,t), 0,h);
        M_stencil(k,l) = int(bk*bl, 0,h);
            
    end
    
end

% now build matrix from stencil
temp = ones(N,1);
K = spdiags([K_stencil(1,2)*temp, (K_stencil(1,1)+K_stencil(2,2))*temp, K_stencil(2,1)*temp], -1:1, N, N);
G = spdiags([G_stencil(1,2)*temp, (G_stencil(1,1)+G_stencil(2,2))*temp, G_stencil(2,1)*temp], -1:1, N, N);
M = spdiags([M_stencil(1,2)*temp, (M_stencil(1,1)+M_stencil(2,2))*temp, M_stencil(2,1)*temp], -1:1, N, N);

% first and last entry are different
K(1,1) = K_stencil(1,1); K(end,end) = K_stencil(2,2);
G(1,1) = G_stencil(1,1); G(end,end) = G_stencil(2,2);
M(1,1) = M_stencil(1,1); M(end,end) = M_stencil(2,2);

end

