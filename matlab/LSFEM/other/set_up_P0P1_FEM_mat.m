%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% build P0, P1 combination %%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [G_p0p1, M_p0p1, M_p0p0] = set_up_P0P1_FEM_mat(N)
h = 1/(N-1);

temp = ones(N,1);
% mass and gradient matrices for basis functions:
% p1 = sum_j \phi_j(t), phi_j linear, hat fct
% p0 = 1 on (t_n, t_n+1) 

% M_p0p1: ii = int(1-1/h*t, 0,h) ii+1 = int(1/h*t, 0,h)
%M_p0p1 = spdiags(0.5*h.*[temp, temp], 0:1, N, N);
M_p0p1 = spdiags(0.5*h.*[temp, temp], 0:1, N-1, N);
%M_p0p1(end, end) = 0;


% M_p0p0: ii = int(1*1, 0,h), zero otherwise
M_p0p0 = spdiags(h*temp, 0, N-1, N-1);
%M_p0p0(1,1) = h; M_p0p0(end, end) = h;

% G_p0p1: ii = int(-1/h, 0,h) ii+1 = int(1/h, 0,h)
G_p0p1 = spdiags([-temp, temp],0:1, N-1, N);
%G_p0p1 = sparse(zeros(N)); G_p0p1(1,1) = -1; G_p0p1(end, end) = 1;



end
