%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% BUILD STANDARD AMG %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% will have to set up problem on FG 
% pre- & post - smoothing
% compute sets of strong dependence
% set up F- and C-nodes using above
% interpolation op: compute weight matrix (w)_ij


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% set up other AMG parameters %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function u = two_level_AMG(A, rhs, u_init, smoother, eps, max_iter)
tol = 10^-6;
max_iter_sm = 3;
% strong dependency matrix, S(i,j)=1 if i str. depends on j
% zero otherwise, C(i)=1 if coarse GP, F(i)=1 if fine GP
S = set_up_str_dep_AMG(A, eps);
[C, F] = set_up_CoarseGr_AMG(S);
save('CoarseGridPts.mat', 'C');

N = length(C);
N_c = sum(C);

% weight matrix W for interpolation
W = set_up_W_AMG(A, C, S);
%restriction will be transpose of interpolation: A_c = I^T*A*I;

%build interpolation operator
I = zeros(N, N_c);
Id = eye(N_c);
count = 1;

for row=1:N
    if(C(row) == 1)
        I(row, :) =  Id(count, :);
        count = count + 1;
    else 
        temp = W(row,:).*S(row,:);
        I(row, :) = temp(C==1);
    end
end

R = I';
A_c = R*A*I;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOOP 
u = u_init;
for k=1:max_iter
    % pre-smoothing
      fprintf('iter : %d\n', k);
      fprintf('energy before pre-smoothing : %d\n', abs(u'*A*u-rhs'*u));
      u = smoothing(A, rhs, u, max_iter_sm, smoother);
      fprintf('energy after pre-smoothing : %d\n', abs(u'*A*u-rhs'*u));
% 
     %u = A \ rhs;
      
    % compute + restrict residual 
     r = rhs - A*u;
     r_c = R*r; 
    
    % solve coarse grid problem directly
     e_c = A_c \ r_c;

    % interpolate back, update u
     e = I*e_c;
     u = u + e;
    
    % post-smoothing
     fprintf('energy before post-smoothing : %d\n', abs(u'*A*u-rhs'*u));
     u = smoothing(A, rhs, u, max_iter_sm, smoother);
     fprintf('energy after post-smoothing : %d\n', abs(u'*A*u-rhs'*u));

    if(norm(A*u-rhs) < tol)
        fprintf('converged after %d iterations \n', k);
        break;
    end
end


end





