function[L_h, R, I, u0, rhs] = set_up_mg(A, f, J, bdy_ind)

% set up Interpolation I/Restriction R matrices of all sizes
% we have J restriction matrices
% first one of size R: (N/2, N), then (N/4, N/2), ... until
% (N/(2^(J-1)), N/(2^J))
% Interpolation matrix
N = length(f);

I = {};
R = {};

% store everything in cell arrays of different sizes, to later perform
% udpdates

L_h = {};
u0 = {};
rhs = {};


% set up matrices for grid transfer 
% only choose every other entry
    for j=1:J+1
    sub_size = floor(N/(2^(j-1)));
    temp = ones(sub_size,1);
    mat_temp = spdiags([temp, 2*temp, temp], -1:1, sub_size, sub_size);
    
   if(j == 1)
    I{j} = (1/2*mat_temp(2:2:end, :))';
    I{j}(bdy_ind, :) = 0;
    %R{j} = 1/4*mat_temp(2:2:end, :);
    R{j} = I{j}';
   end
   
    I{j} = (1/2*mat_temp(2:2:end, :))';
    %R{j} = 1/4*mat_temp(2:2:end, :);
    R{j} = I{j}';
    
    % to store intermediate results
    u0{j} = zeros(sub_size,1);  
    
    if j==1
        L_h{1} = A;
        rhs{1} = f;
    else
        
    L_h{j} = R{j-1}*L_h{j-1}*I{j-1};
    rhs{j} = zeros(sub_size, 1);
    end
    
    end

    end