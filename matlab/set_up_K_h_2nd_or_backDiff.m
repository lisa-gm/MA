function K_h = set_up_K_h_backDiff(Nt, Nx, c, delta_t, delta_x)
%%% I assume that we are on [0,1], [0,1]^2 

    % only consider inner grid points 
    % otherwise rows of 1's in there in between
    pts = Nx;
    total_pts = Nt*Nx;   
%% what if I use a different difference scheme that only uses 
    % backward differences
    % => won't be symmetric anymore
    
    % 2 -5 4 -1
    
    temp = ones(total_pts,1);

    B = spdiags([-c/delta_x^2*temp, (2*c/delta_x^2+2/delta_t^2)*temp, -c/delta_x^2*temp], -1:1, pts, pts);
    B_list = {};
    
    for i=1:Nt
        B_list{i} = B;
    end
    
    K_h_blk = sparse(blkdiag(B_list{1:end}));
    
    K_h_sides =  sparse(spdiags([-5/delta_t^2*temp, 4/delta_t^2*temp, -1/delta_t^2*temp], [-(pts),-(2*pts), -3*pts], total_pts, total_pts));
    K_h = (K_h_blk + K_h_sides);
    
    % the first 2 blocks cannot change, we have the initial conditions,
    % hence replace first row block by identity and then zeros
    
    K_h(1:3*pts,1:3*pts) = eye(3*pts);
    K_h(1:3*pts,3*pts+1:end) = zeros(3*pts, total_pts - 3*pts);
    
    % we also have fixed boundary conditions for first and last entry
    % of every block in vector, make rows with entry 1
    
    K_h(1:pts:end, :) = 0;
    K_h(pts:pts:end, :) = 0; 
    
    K_h(1:pts:end, 1:pts:end) = eye(Nt);
    K_h(pts:pts:end, pts:pts:end) = eye(Nt); 
    
    %spy(K_h)
    
end