function K_global = set_up_FEM_2D_stiffness(Nx, Nt)
hx = 1/(Nx-1);
ht = 1/(Nt-1);

offset = Nx;

%% time operators
%syms t x

%% basis function on reference element q = 1 

grad_l1l1 =  @(x,t) 1/hx^2*(1-t./ht).^2+ (1-x./hx).^2.*(-1/ht)^2;
grad_l1l2 =  @(x,t) (-1/hx).*(1-t./ht).*  1/hx .*(1-t./ht) + x.*x./hx^2 .*(-1/ht)^2;
grad_l1l3 =  @(x,t) (-1/hx).*(1-t./ht).* (-1/hx).* t./ht + (1-x./hx).*(-1/ht).* (1-x./hx).* 1/ht;
grad_l1l4 =  @(x,t) (-1/hx).*(1-t./ht).*  1/hx .* t./ht + (1-x./hx).*(-1/ht).*  x./hx .* 1/ht;

grad_l2l2 =  @(x,t)  1/hx^2 .*(1-t./ht).^2  +  (x./hx).^2 .*(1/ht)^2;
grad_l2l3 =  @(x,t)  1/hx .*(1-t./ht)  .*  (-1/hx).* t./ht + x./hx .*(-1/ht)  .*  (1-x./hx).* 1/ht;
grad_l2l4 =  @(x,t)  1/hx .*(1-t./ht)  .*  1/hx .* t./ht + x./hx .*(-1/ht)  .*  x./hx .* 1/ht;

grad_l3l3 =  @(x,t)  (-1/hx)^2.* (t./ht).^2  + (1-x./hx).^2.* 1/ht^2;
grad_l3l4 =  @(x,t)  (-1/hx).* t./ht  .*  1/hx .* t./ht +  (1-x./hx).* 1/ht  .*  x./hx .* 1/ht;

grad_l4l4 =  @(x,t)   1/hx^2 .* (t./ht).^2    +   (x./hx).^2 .* 1/ht^2;

K_local = zeros(4);
K_global = zeros(Nx*Nt);

%%%%%%%%%%%%%%%%%%%% LOCAL MASS MATRIX %%%%%%%%%%%%%%%%%%%%%%%%
K_local(1,1) = integral2(grad_l1l1, 0, hx, 0, ht);
K_local(1,2) = integral2(grad_l1l2, 0, hx, 0, ht);
K_local(1,3) = integral2(grad_l1l3, 0, hx, 0, ht);
K_local(1,4) = integral2(grad_l1l4, 0, hx, 0, ht);

K_local(2,2) = integral2(grad_l2l2, 0, hx, 0, ht);
K_local(2,3) = integral2(grad_l2l3, 0, hx, 0, ht);
K_local(2,4) = integral2(grad_l2l4, 0, hx, 0, ht);

K_local(3,3) = integral2(grad_l3l3, 0, hx, 0, ht);
K_local(3,4) = integral2(grad_l3l4, 0, hx, 0, ht);

K_local(4,4) = integral2(grad_l4l4, 0, hx, 0, ht);

K_local(2,1) = K_local(1,2); K_local(3,1) = K_local(1,3); K_local(4,1) = K_local(1,4);
K_local(3,2) = K_local(2,3); K_local(4,2) = K_local(2,4);
K_local(4,3) = K_local(3,4);

K_local 

for J=1:Nt-1
    for I=1:Nx-1
        local2Global=[ (J-1)*offset+I (J-1)*offset+I+1 (J)*offset+I (J)*offset+I+1];
        
        for i=1:4
            for j=1:4
                K_global(local2Global(i),local2Global(j))=K_global(local2Global(i),local2Global(j))+K_local(i,j);
            end
        end
        
        
    end
end

end