function x = GaussRadau_RightPoint(N, a, b)

% Nodes are roots of Jacobi polynomial J^(1,0)_(N-1)
% see http://de.wikipedia.org/wiki/Jacobi-Polynom

alpha = 1.0;
beta  = 0.0;

diag    = zeros(N-1,1);
subdiag = zeros(N-2,1);

diag(1)    = (beta - alpha)/(2+alpha+beta);

for jj=1:N-2
    diag(jj+1)   = (beta-alpha)*(alpha+beta)/(2*jj + 2 + alpha + beta)/(2*jj+alpha+beta);
    subdiag(jj)  = sqrt( 4*jj*(jj+alpha)*(jj+beta)*(jj+alpha+beta) )/sqrt( (2*jj-1+alpha+beta)*(2*jj+alpha+beta)^2*(2*jj+1+alpha+beta)); 
end

subdiag1 = [subdiag; 0];
subdiag2 = [0; subdiag];

Mat = full(spdiags([subdiag1 diag subdiag2], -1:1, N-1, N-1));

x = sort(eig(Mat));

if alpha==1
    x = [x ; 1.0];
elseif beta==1
    x = [-1.0; x];
end

% Linear map from[-1,1] to [a,b] (from GaussLegendre)
x=(a*(1-x)+b*(1+x))/2;

x = x.';