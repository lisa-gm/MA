function L_h = set_up_L_h_FD(N)
%%% I assume that we are on [0,1], [0,1]^2 or [0,1]^3 respectively

% do two versions of only inner system or with boundary as well
pts = N;
h = 1/(N-1);
K = 1;

temp = ones(pts,1);
L_h = K/h^2 * spdiags([-temp, 2*temp, -temp], -1:1, pts, pts);

% zero boundary conditions PUT BACK IN
L_h(1,:) = [1, zeros(1,N-1)];
L_h(end,:) = [zeros(1,N-1), 1];
% also eliminate entry in L_h(2,1) and L_h(N, N+1) to keep symmetry
% hence we have to update f(2) and f(N) later!
%L_h(2,1) = 0;
%L_h(N,N+1) = 0;
end 