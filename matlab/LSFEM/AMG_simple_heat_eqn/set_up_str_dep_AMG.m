%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% compute sets of strong dependence %%%%%%%

% S(i,j) = 1 if i strongly depends on j

function S = set_up_str_dep_AMG(A, eps)
tic;
N = size(A,1);
S = zeros(N);
A_no_diag = A - diag(diag(A));

% if row only one entry, all other zero, assume value predefined
% set whole row to zero, does not need to be in coarse grid!

for i=1:N
    max_a = max(abs(A_no_diag(i,:)));
    for j=1:N
        if(and(-A(i,j) >= eps*max_a, max_a > 0))
            S(i,j) = 1;
        end
    end
end

t = toc;
% fprintf('set up strong dep matrix finished, time: %d sec\n', t);

end