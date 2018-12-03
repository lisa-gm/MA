function K_h = set_up_K_h_backDiff(Nt, Nx, c, delta_t, delta_x)
%%% points of the same time will be next to each other
    % while points of different time but same space will have distance
    % multiple of Nx

    % only consider inner grid points 
    % otherwise rows of 1's in there in between
    pts = Nx;
    total_pts = Nt*Nx;   

%% what if I use a different difference scheme that only uses 
    % backward differences
    % => won't be symmetric anymore
    
    % entries in the main diagonal blocks
    temp = ones(total_pts,1);
    B = spdiags([-c/delta_x^2*temp, (2*c/delta_x^2+1/delta_t^2)*temp, -c/delta_x^2*temp], -1:1, pts, pts);
    B_list = {};
    
    for i=1:Nt
        B_list{i} = B;
    end
    
    K_h_blk = blkdiag(B_list{1:end});
    
    % set up the two lower diagonal blocks
    K_h_sides =  spdiags([-2/delta_t^2*temp, 1/delta_t^2*temp], [-(pts),-(2*pts)], total_pts, total_pts);
    K_h = (K_h_blk + K_h_sides);
    
    % the first 2 blocks cannot change, we have the initial conditions,
    % hence replace first row block by identity and then zeros
    
    K_h(1:2*pts,1:2*pts) = eye(2*pts);
    K_h(1:2*pts,2*pts+1:end) = zeros(2*pts, total_pts - 2*pts);
    
    % we also have fixed boundary conditions for first and last entry
    % of every block in vector, make rows with entry 1
    
    K_h(1:pts:end, :) = 0;
    K_h(pts:pts:end, :) = 0; 
    
    K_h(1:pts:end, 1:pts:end) = eye(Nt);
    K_h(pts:pts:end, pts:pts:end) = eye(Nt); 
    
    spy(K_h)
    
end