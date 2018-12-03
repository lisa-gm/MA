%
% Tests the functionality of the inter_poly class.
%
clearvars -except test_all_dir test_all_ll test_all_classname

addpath('../src');

f = @(x) sin(pi*x);
% For too high order, the condition of the problem will cause the test for
% 1e-12 satisfaction of the interpolation condiiton to fail
nmax = 15;
err  = zeros(1,nmax);

for nn=1:nmax
    xi = linspace(0,0.5,nn);
    fi = f(xi);
    pol = interp_poly(xi,fi);
    
    % First make sure that the polynomial really interpolates
    fi_pol = pol.evaluate(xi);
    assert( norm(fi - fi_pol, inf) < 1e-12, 'Polynomial does not interpolate');
       
end

% Test that the interpolation polynomial reproduces monomials. For higher
% order the computations suffer from a noticeable accumulation of round of
% error, so the bound of 1e-10 set below won't be reached anymore.
for nn=1:9
    x_test = linspace(0, 0.5, 1000);
    for pp=1:(nn-1)
        f = @(x) x.^pp;
        xi    = linspace(0,5,nn);
        fxi   = f(xi);
        % Because f naturally interpolates (xi, fxi) and the interpolation
        % polynomial is unique, the resulting interp_poly should be
        % identical to f
        pol   = interp_poly(xi, fxi);
        
        % Verify that f and interp_poly yield (approximately) identical
        % results on a number of test points (x_test)
        y_pol = pol.evaluate(x_test);
        assert( norm(y_pol - f(x_test), inf) < 1e-10, 'Polynomial failed to reproduce monomial with order %2d', pp);
    end
end
    
fprintf('[0] - INTERP_POLY - All tests successful \n');