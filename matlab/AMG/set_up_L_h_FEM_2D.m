%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% 5 POINT STENCIL FEM 2D %%%%%%%%%%%%%%%%%%%%%

% please note that for hom. Neumann bdy conditions: "no bdy"
% grid becomes "1 larger" 

function L_h = set_up_L_h_FEM_2D(Nx, Nt, bdy_cond)
%%% I assume that we are on [0,1], [0,1]^2 or [0,1]^3 respectively

pts = Nx;
h = 1/(pts-1);

L_h = zeros(Nt*pts);
temp = ones(pts*Nt,1);
L_h_blk = spdiags([-temp, 4*temp, -temp], -1:1, pts, pts);

if(strcmp(bdy_cond, 'Dirichlet'))
% zero boundary conditions PUT BACK IN
L_h_blk(1,:) = [1, zeros(1,Nx-1)];
L_h_blk(end,:) = [zeros(1,Nx-1), 1];
% also eliminate entry in L_h(2,1) and L_h(N, N+1) to keep symmetry
% hence we have to update f(2) and f(N) later!
%L_h(2,1) = 0;
%L_h(N,N+1) = 0;

% for off-diagonal part
bdy_pts_w = 1:Nx:Nx*Nt;
bdy_pts_e = Nx:Nx:Nx*Nt;
temp(bdy_pts_w) = 0;
temp(bdy_pts_e) = 0;
end 


for i=1:Nt
    L_h((i-1)*Nx+1:i*Nx, (i-1)*Nx+1:i*Nx) = L_h_blk;
end

L_h_sides = spdiags([-temp, -temp], [-pts, pts], pts*Nt, pts*Nt);
L_h = 1/h^2*(sparse(L_h) + L_h_sides);

end 