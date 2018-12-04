%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% setting up COARSE GRID POINTS %%%%%%%%%%%%%%%%

% compute \lambda parameter: 
% lambda_i = |S_i^T \intersect U| + 2|S_i^T \interset F|
% then select coarse grid points

function [C, F] = set_up_CoarseGr_AMG(S)
tic;
N = size(S, 1);

U = ones(1, N);     % first all gridpts are undecided
F = zeros(1, N);   % no fine gridpts
C = zeros(1, N);

for i=1:N
    if(sum(S(i,:)) == 0)
        F(i) = 1;
        U(i) = 0;
    end
end


counter = 0;
while(sum(U) > 0)
    lam = zeros(N, 1);
    ind_one = find(U);
    
    for it=1:sum(U)
        lam(ind_one(it)) = U*S(:,ind_one(it)) + 2*(F*(S(:,ind_one(it)))); 
    end
    
    [max_lam, max_lam_ind] = max(lam);
    
    if(max_lam == 0)
        max_lam_ind = ind_one(1);
    end
    
    % max_lam_ind new C coarse grid pt, take off U
    C(max_lam_ind) = 1;
    U(max_lam_ind) = 0;
    % now find all those strongly connected to max_lam_ind in U
    F = F + (S(:,max_lam_ind)'.*U);
    U(S(:,max_lam_ind) == 1) = 0;
    
    counter = counter +1; 
    if(counter > N)
        fprintf('something wrong in set up coarse grid AMG\n');
        
        break;
    end
end
t = toc;
% fprintf('set up caorse grid finished, time: %d sec\n', t);
end