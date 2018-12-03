%% setting up the DG matrices for time stepping

function Lh_loc = set_up_Lh_FEM_2D_local_basis(hx)
%% time operators
syms x y

%% basis function on reference element q = 1 
l1(x,y) = x/hx;
l2(x,y) = y/hx;
l3(x,y) = 1-x/hx-y/hx;
basis = {l1,l2,l3};

Lh_loc = zeros(3);


for k = 1:3
    for l = 1:3
        
        bk = basis{k};
        bl = basis{l};
        
        ymax = @(x) 1 - x;
        % q = integral2(fun,0,1,0,ymax)
        
        Lh_loc(k,l) = integral2(diff(bk,x)*diff(bl,x)+diff(bk,y)*diff(bl,y),0,hx,0,ymax);
        % K_t(k,l) = - int(diff(bk,t)*bl,0,dt) + bk(dt)*bl(dt);              
               
    end
end