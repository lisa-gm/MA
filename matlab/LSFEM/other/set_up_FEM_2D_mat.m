function M_local = set_up_FEM_2D_mat(Nx, Nt)
hx = 1/(Nx-1);
ht = 1/(Nt-1);

%% time operators
%syms t x

%% basis function on reference element q = 1 
l1 = @(x,t) (1-x./hx).*(1-t./ht);
l2 = @(x,t)   x./hx .*(1-t./ht);
l3 = @(x,t) (1-x./hx).* t./ht;
l4 = @(x,t)    x./hx .* t./ht;
basis = {l1,l2,l3,l4};

l1l1 =  @(x,t) (1-x./hx).*(1-t./ht).* (1-x./hx).*(1-t./ht);
l1l2 =  @(x,t) (1-x./hx).*(1-t./ht).*  x./hx .*(1-t./ht);
l1l3 =  @(x,t) (1-x./hx).*(1-t./ht).* (1-x./hx).* t./ht;
l1l4 =  @(x,t) (1-x./hx).*(1-t./ht).*  x./hx .* t./ht;

l2l2 =  @(x,t)  x./hx .*(1-t./ht)  .*  x./hx .*(1-t./ht);
l2l3 =  @(x,t)  x./hx .*(1-t./ht)  .*  (1-x./hx).* t./ht;
l2l4 =  @(x,t)  x./hx .*(1-t./ht)  .*  x./hx .* t./ht;

l3l3 =  @(x,t)  (1-x./hx).* t./ht  .* (1-x./hx).* t./ht;
l3l4 =  @(x,t)  (1-x./hx).* t./ht  .*  x./hx .* t./ht;

l4l4 =  @(x,t)   x./hx .* t./ht    .*   x./hx .* t./ht;

% @(x,t)

%G_t_stencil = zeros(4);
%G_h_stencil = zeros(4);
K_local = zeros(4);
%K_sp_stencil = zeros(4);
M_local = zeros(4);

%%%%%%%%%%%%%%%%%%%% LOCAL MASS MATRIX %%%%%%%%%%%%%%%%%%%%%%%%
M_local(1,1) = integral2(l1l1, 0, hx, 0, ht);
M_local(1,2) = integral2(l1l2, 0, hx, 0, ht);
M_local(1,3) = integral2(l1l3, 0, hx, 0, ht);
M_local(1,4) = integral2(l1l4, 0, hx, 0, ht);

M_local(2,2) = integral2(l2l2, 0, hx, 0, ht);
M_local(2,3) = integral2(l2l3, 0, hx, 0, ht);
M_local(2,4) = integral2(l2l4, 0, hx, 0, ht);

M_local(3,3) = integral2(l3l3, 0, hx, 0, ht);
M_local(3,4) = integral2(l3l4, 0, hx, 0, ht);

M_local(4,4) = integral2(l4l4, 0, hx, 0, ht);

M_local(2,1) = M_local(1,2); M_local(3,1) = M_local(1,3); M_local(4,1) = M_local(1,4);
M_local(3,2) = M_local(2,3); M_local(4,2) = M_local(2,4);
M_local(4,3) = M_local(3,4);

%%%%%%%%%%%%%%%%% LOCAL STIFFNESS MATRIX %%%%%%%%%%%%%%%%

% gl1gl1 =  @(x,t) (1/hx).*(1-t./ht).* (1-x./hx).*(1-t./ht);
% gl1gl2 =  @(x,t) (1-x./hx).*(1-t./ht).*  x./hx .*(1-t./ht);
% gl1gl3 =  @(x,t) (1-x./hx).*(1-t./ht).* (1-x./hx).* t./ht;
% gl1gl4 =  @(x,t) (1-x./hx).*(1-t./ht).*  x./hx .* t./ht;
% 
% gl2gl2 =  @(x,t)  x./hx .*(1-t./ht)  .*  x./hx .*(1-t./ht);
% gl2gl3 =  @(x,t)  x./hx .*(1-t./ht)  .*  (1-x./hx).* t./ht;
% gl2gl4 =  @(x,t)  x./hx .*(1-t./ht)  .*  x./hx .* t./ht;
% 
% gl3gl3 =  @(x,t)  (1-x./hx).* t./ht  .* (1-x./hx).* t./ht;
% gl3gl4 =  @(x,t)  (1-x./hx).* t./ht  .*  x./hx .* t./ht;
% 
% gl4gl4 =  @(x,t)   x./hx .* t./ht    .*   x./hx .* t./ht;

%%%%%%%%%%% ASSEMBLE GLOBAL MASS MATRIX %%%%%%%%%%%%%%%%%%%%%

M = spdiags([M_stencil(2,1)*temp, (M_stencil(1,1)+M_stencil(2,2))*temp, M_stencil(1,2)*temp], -1:1, Nx*Nt, Nx*Nt);


% temp = ones(N,1);
% K = spdiags([K_stencil(2,1)*temp, (K_stencil(1,1)+K_stencil(2,2))*temp, K_stencil(1,2)*temp], -1:1, N, N);
% G = spdiags([G_stencil(2,1)*temp, (G_stencil(1,1)+G_stencil(2,2))*temp, G_stencil(1,2)*temp], -1:1, N, N);
% M = spdiags([M_stencil(2,1)*temp, (M_stencil(1,1)+M_stencil(2,2))*temp, M_stencil(1,2)*temp], -1:1, N, N);
% 
% % first and last entry are different
% K(1,1) = K_stencil(1,1); K(end,end) = K_stencil(2,2);
% G(1,1) = G_stencil(1,1); G(end,end) = G_stencil(2,2);
% M(1,1) = M_stencil(1,1); M(end,end) = M_stencil(2,2);



end