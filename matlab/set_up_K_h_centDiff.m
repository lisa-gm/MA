%%%%%%%%%%% setting up the finite difference stencil %%%%%%%%%%%%%

% for now implementation only works for dim = 1 !! 

function K_h = set_up_K_h_centDiff(Nt, Nx, p, delta_x, delta_t)
    
    total_pts = (Nx+1)*(Nt+1);
    pts = Nx+1;
    
    % time to set up 5-point stencil 
    temp = ones(total_pts,1);

    B = spdiags([-p/delta_x^2*temp, (2*p/delta_x^2-2/delta_t^2)*temp, -p/delta_x^2*temp], -1:1, Nx+1, Nx+1);
    B_list = {};
    
    for i=1:Nt+1
        B_list{i} = B;
    end
    
    K_h_blk = blkdiag(B_list{1:end});
    
    K_h_sides =  spdiags([1/delta_t^2*temp, 1/delta_t^2*temp], [-(Nx-1),(Nx-1)], total_pts, total_pts);
    K_h = (K_h_blk + K_h_sides);
    
    % the first block cannot change, we have the initial conditions,
    % hence replace first row block by identity and then zeros
    
    K_h(1:Nx+1,1:Nx+1) = eye(Nx+1);
    K_h(1:Nx+1,Nx+2:end) = zeros(Nx+1, total_pts - Nx-1);
    
    K_h(end-Nx:end,end-Nx:end) = eye(Nx+1);
    K_h(end-Nx:end,1:end-Nx-1) = zeros(Nx+1, total_pts - Nx-1);

    
    % we also have fixed boundary conditions for first and last entry
    % of every block in vector, make rows with entry 1
    
    K_h(1:pts:end, :) = 0;
    K_h(pts:pts:end, :) = 0; 
    
    K_h(1:pts:end, 1:pts:end) = eye(Nt+1);
    K_h(pts:pts:end, pts:pts:end) = eye(Nt+1);
    

end
