%%%%%%%%%%% setting up the finite element stencil %%%%%%%%%%%%%

% for now implementation only works for dim = 1 !! 

function K_h = set_up_K_h(Nt, Nx, p, dt, dx)
%%% I assume that we are on [0,1], [0,1]^2 
    %s = p*(dt/dx)^2;
    % only consider inner grid points 
    % otherwise rows of 1's in there in between
    pts = Nx;
    total_pts = Nt*Nx;

    temp = ones(total_pts,1);

    B = spdiags([-p/dx^2*temp, (2*p/dx^2-2/dt^2)*temp, -p/dx^2*temp], -1:1, pts, pts);
    %B = spdiags([-temp, 2-2*s*temp, -temp], -1:1, pts, pts);
    B_list = {};
    
    for i=1:Nt
        B_list{i} = B;
    end
    
    K_h_blk = blkdiag(B_list{1:end});
    
    K_h_sides =  spdiags([1/dt^2*temp, 1/dt^2*temp], [-Nx,Nx], total_pts, total_pts);
    %K_h_sides =  spdiags([s*temp, s*temp], [-N,N], total_pts, total_pts);

    K_h = (K_h_blk + K_h_sides);
        
    % the first block cannot change, we have the initial conditions,
    % hence replace first row block by identity and then zeros
    
    K_h(1:pts,1:pts) = eye(pts);
    K_h(1:pts,pts+1:end) = zeros(pts, total_pts - pts);
    
    K_h(total_pts-pts+1:end,total_pts-pts+1:end) = eye(pts);
    K_h(total_pts-pts+1:end,1:total_pts-pts) = zeros(pts, total_pts - pts);
    
    % we also have fixed boundary conditions for first and last entry
    % of every block in vector, make rows with entry 1
    
    K_h(1:pts:end, :) = 0;
    K_h(pts:pts:end, :) = 0; 
    
    K_h(1:pts:end, 1:pts:end) = eye(Nt);
    K_h(pts:pts:end, pts:pts:end) = eye(Nt);
    
    %spy(K_h)
    
    %% what if I use a different difference scheme that only uses 
    % backward differences
    % => won't be symmetric anymore
    
    % 2 -5 4 -1
    
%     temp = ones(total_pts,1);
% 
%     B = spdiags([-p/delta_x^2*temp, (2*p/delta_x^2+2/delta_t^2)*temp, -p/delta_x^2*temp], -1:1, pts, pts);
%     B_list = {};
%     
%     for i=1:pts
%         B_list{i} = B;
%     end
%     
%     K_h_blk = blkdiag(B_list{1:end});
%     
%     K_h_sides =  spdiags([-5/delta_t^2*temp, 4/delta_t^2*temp, -1/delta_t^2*temp], [-(N-1),-2*N+1, -3*N+1], total_pts, total_pts);
%     K_h = (K_h_blk + K_h_sides);
%     
%     % the first block cannot change, we have the initial conditions,
%     % hence replace first row block by identity and then zeros
%     
%     K_h(1:pts,1:pts) = eye(pts);
%     K_h(1:pts,pts+1:end) = zeros(pts, total_pts - pts);
%     
%     % we also have fixed boundary conditions for first and last entry
%     % of every block in vector, make rows with entry 1
%     
%     K_h(1:pts:end, :) = 0;
%     K_h(pts:pts:end, :) = 0; 
%     
%     K_h(1:pts:end, 1:pts:end) = eye(pts);
%     K_h(pts:pts:end, pts:pts:end) = eye(pts); 
    

end
