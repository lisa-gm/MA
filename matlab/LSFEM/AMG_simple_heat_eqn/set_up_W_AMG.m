%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% SET UP INTERPOLATION OPERATOR AMG %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% we have 
% -- C: 0 if FG, 1 if CG point
% -- F: 1 if FG, 0 if CG point
% -- now: Cs (i,j) = 1 if j is a COARSE GRID POINT & i strongly 
% depends on j, compute by componentwise mult. of 
%  S(i, :).*C
% -- Ds(i,j) = 1 if j is a FINE GP and i strongly depends on j
% compute: S(i,:).*F
% -- Dw = ones(N) - Cs - Ds

% ith comp of (Ih_2h e)_i = e_i if i \in C or
% sum_j \in C_i w_ij*e_j if i \in F

% compute w_ij now

function W = set_up_W_AMG(A, C, S)
tic;
N = size(A, 1);
W = zeros(N);

% fine grid points
F = ones(1,N) - C;

Cs = S.*C;
Ds = S.*F;
Dw = ones(N) - Cs - Ds;

AtimesC = A*Cs'; 

for i=1:N
    for j=1:N        
    temp_par = (Ds(i,:).*A(i,:).*A(:,j)')./AtimesC(:,i)';
    temp_par(AtimesC(:,i)==0) = 0;
    temp_sum = sum(temp_par); 

    temp_down = dot(Dw(i,:),A(i,:));
    W(i,j) = (A(i,j) + temp_sum)/(A(i,i) + temp_down);
    end
end
t = toc;

row_sum_W = sum(W,2);
% fprintf('set up weight matrix finished, time: %d sec\n', t);
end 