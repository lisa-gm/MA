%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% FEM MASS MATRIX 2D %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% local matrix looks as follows: h^2/24* [2, 1, 1; 1, 2, 1; 1, 1, 2]
% double check h^2
% 3 different blocks, diag, upper, lower

function M_h = set_up_L_h_FEM_2D(Nx, Nt, bdy_cond)

hx = 1/(Nx-1);
temp = ones(Nx, 1);

M_h_diag = spdiags([temp, 2*temp, temp], -1:1, Nx, Nx);
M_h_low = spdiags([temp, temp, 2*temp], -1:1, Nx, Nx);
M_h_up = spdiags([2*temp, temp, temp], -1:1, Nx, Nx);

if(strcmp(bdy_cond, 'Dirichlet'))
    M_h_diag(1,:) = [1, zeros(1, Nx-1)];
    M_h_diag(end,:) = [zeros(1, Nx-1), 1];
    
    M_h_low(1, :) = 0;
    M_h_low(end, :) = 0;
    
    M_h_up(1, :) = 0;
    M_h_up(end, :) = 0;

end

M_h = zeros(Nx*Nt);

for i=1:Nt
    M_h((i-1)*Nx+1:i*Nx, (i-1)*Nx+1:i*Nx) = M_h_diag;
    
    if(i>1)
        M_h((i-2)*Nx+1:(i-1)*Nx, (i-1)*Nx+1:i*Nx) = M_h_up; end
    
    if(i<Nt)
        M_h(i*Nx+1:(i+1)*Nx, (i-1)*Nx+1:i*Nx) = M_h_low; end
    
end

M_h = hx^2/24*sparse(M_h);

end