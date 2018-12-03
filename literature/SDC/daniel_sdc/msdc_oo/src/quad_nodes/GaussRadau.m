function [x,w] = GaussRadau(N,a,b)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% lgrnodes.m
%
% Computes the Legendre-Gauss-Radau nodes, weights and the LGR Vandermonde
% matrix. The LGR nodes are the zeros of P_N(x)+P_{N+1}(x).
%
% References on LGR nodes and weights:
%   C. Canuto, M. Y. Hussaini, A. Quarteroni, T. A. Tang, "Spectral Methods
%   in Fluid Dynamics," Section 2.3. Springer-Verlag 1987
%
%   F. B. Hildebrand , "Introduction to Numerical Analysis," Section 8.11
%   Dover 1987
%
% Written by Greg von Winckel - 05/02/2004
% Contact: gregvw@chtm.unm.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% subtract 1 from N so that length(x)==nr_nodes
N=N-1;

% Truncation + 1
N1=N+1;

% Use Chebyshev-Gauss-Radau nodes as initial guess for LGR nodes
x=-cos(2*pi*(0:N)/(2*N+1))';

% The Legendre Vandermonde Matrix
P=zeros(N1,N1+1);

% Compute P_(N) using the recursion relation
% Compute its first and second derivatives and
% update x using the Newton-Raphson method.

xold=2;

% Free abscissae
free=2:N1;

% Solve F(x) = P_N(x) + P_{N+1}(X) = 0 using Newton
while max(abs(x-xold))>eps
    
    xold=x;
    
    % P_N(-1) = (-1)^N for all Legendre polynomials P_N
    P(1,:)=(-1).^(0:N1);
    
    % P(:,1) = P_0(x) = 1
    P(free,1)=1;    
    
    % P(:,2) = P_1(x) = x
    P(free,2)=x(free);
    
    % Compute P(:,k) for 3,...,N1+1 using the Legendre polynomial recurrence
    for k=2:N1
        % Recurrence: P_N(x) = [ (2N-1)*x*P_(N-1)(x) - (N-1)*P_(N-2)(x) ]/N
        P(free,k+1)=( (2*k-1)*x(free).*P(free,k)-(k-1)*P(free,k-1) )/k;
    end
    
    % Derivative: P_N'(x) = N*( P_(N-1)(x) - x*P_N(x) )/(1-x^2)
    
    % F(x)  = P(free,N1) + P(free,N1+1)
    % DF(x) = P'(free,N1) + P'(free,N1+1) = 

    % xnew = xold - (1/DF(x))*F(x)
    x(free)=xold(free)-((1-xold(free))/N1).*(P(free,N1)+P(free,N1+1))./(P(free,N1)-P(free,N1+1));
end


% The Legendre-Gauss-Radau Vandermonde
P=P(1:N1,1:N1);

% Compute the weights
w=zeros(N1,1);
w(1)=2/N1^2;
w(free)=(1-x(free))./(N1*P(free,N1)).^2;

% DIRTY HACK

% Linear map from[-1,1] to [a,b] (from GaussLegendre)
x=(a*(1-x)+b*(1+x))/2;

% Scale weights according to derivative of transformation
w = ((b-a)/2)*w;

x = x.';
w = w.';

end
