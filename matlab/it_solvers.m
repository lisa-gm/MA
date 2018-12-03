N = 100; 
% add constant K
K = 1;

h = 1/N;
x = linspace(0,1, N+1);

%%% here we know the real solution: 
u_real = 1/(4*pi^2) * sin(2*pi*x);

L_h = K/h^2 * full(gallery('tridiag',N+1,-1,2,-1));    
L_h(1,:) = [1, zeros(1,N)];
L_h(end,:) = [zeros(1,N), 1];
% also eliminate entry in L_h(2,1) and L_h(N, N+1) to keep symmetry
% hence we have to update f(2) and f(N) later! 
L_h(2,1) = 0;
L_h(N,N+1) = 0;

% and the RHS f
f = transpose(sin(2*pi*x));

% put in boundary conditions in case we change f
f(1) = 0;
f(end) = 0;

% and how the effect the neighbouring values
f(2) = f(2) - K/h^2 * f(1);
f(end-1) = f(end-1) - K/h^2 * f(end);


% system we would like to solve L_h * u = f for u
u0 = zeros(N+1,1);

% call Jacobi Solve & Gauss Seidel Solve
[u_J, error_J] = JacobiSolve(L_h, f, u0, 3000);
[u_GS, error_GS] = GaussSeidelSolve(L_h, f, u0, 3000);

figure
plot(x, u_J, x, u_GS, x, u_real)
legend('u Jacobi', 'u Gauss Seidel', 'u real');
