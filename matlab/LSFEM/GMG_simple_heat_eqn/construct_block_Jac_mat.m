% *********************************************************************** %
% ***************** BLOCK JACOBI SMOOTHER EXTENDED! ********************* %
% *********************************************************************** %

% considering blocks of size block_size for sigma and u for same dof on
% grid

% here we walk in time then space, could also couple time or do patches 
% of space-time

% assemble whole thing 

function P = construct_block_Jac_mat(Nx_pts, Nt_pts, A, block_size_s, block_size_t)
tot_pts = Nx_pts*Nt_pts;

% since block size can vary maybe not divisible
% careful to not couple grid(i, end) and grid(i+1, 1) ... 
k_s = floor(Nx_pts/block_size_s);
last_mult_s = k_s*block_size_s;
rest_size_s = Nx_pts - last_mult_s;

k_t = floor(Nt_pts/block_size_t);
last_mult_t = k_t*block_size_t;
rest_size_t = Nt_pts - last_mult_t;

P = zeros(size(A));

    for j = 1:block_size_t:last_mult_t
        for i=1:block_size_s:last_mult_s
            
            curr_inds = [];
            for l=0:block_size_t-1
                curr_inds = [curr_inds, (j-1+l)*Nx_pts+i: (j-1+l)*Nx_pts+i-1+block_size_s, (j-1+l)*Nx_pts+i+tot_pts: (j-1+l)*Nx_pts+i-1+block_size_s+tot_pts];   
            end
          
            curr_inds = sort(curr_inds);
           
            sub_M =  inv(A(curr_inds, curr_inds));
            P(curr_inds, curr_inds) = sub_M;
        end
               
        % rest size indices never reached ... for specific t, update now
        if(rest_size_s ~= 0)
            curr_inds = [];
            
            for l=0:block_size_t-1
                curr_inds = [curr_inds, (j-1+l)*Nx_pts+[last_mult_s+1:Nx_pts], tot_pts+(j-1+l)*Nx_pts+[last_mult_s+1:Nx_pts]];
            end
            curr_inds = sort(curr_inds);
            
            sub_M =  inv(A(curr_inds, curr_inds));
            P(curr_inds, curr_inds) = sub_M;
        end
          
    end
    
    if(rest_size_t ~= 0)
       
    for i=1:block_size_s:last_mult_s
    curr_inds = [];
        for l=0:rest_size_t-1
            curr_inds = [curr_inds, (last_mult_t+l)*Nx_pts+i: (last_mult_t+l)*Nx_pts+i-1+block_size_s, (last_mult_t+l)*Nx_pts+i+tot_pts: (last_mult_t+l)*Nx_pts+i-1+block_size_s+tot_pts];
        end
    curr_inds = sort(curr_inds);
            
    sub_M =  inv(A(curr_inds, curr_inds));
    P(curr_inds, curr_inds) = sub_M;
    end
              
    curr_inds = [];
            
    for l=0:rest_size_t-1
        curr_inds = [curr_inds, (last_mult_t+l)*Nx_pts+[last_mult_s+1:Nx_pts], tot_pts+(last_mult_t+l)*Nx_pts+[last_mult_s+1:Nx_pts]];
    end
    curr_inds = sort(curr_inds);
            
    sub_M =  inv(A(curr_inds, curr_inds));
    P(curr_inds, curr_inds) = sub_M;
    end
    
end