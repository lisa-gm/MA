%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% FINDING GOOD INITIAL GUESS NEWTON %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function u_init = init_guess_newton(f, df, d2f, u0, max_iter)

F = @(v) u0 + df(v-u0) - f(v);
dF = @(v) diag(d2f(v-u0) - df(v));
for it=1:max_iter
    u0 = u0 - dF(u0) \ F(u0);
end
u_init = u0;
F(u_init);
end