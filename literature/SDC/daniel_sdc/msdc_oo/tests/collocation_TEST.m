% 
% Tests the collocation class
%
clearvars -except test_all_dir test_all_ll test_all_classname

addpath('../src');
addpath('../src/quad_nodes');

max_nodes = [15, 15, 15, 15, 15, 15, 12];
type      = {'gauss-chebyshev', 'gauss-legendre', 'gauss-lobatto', 'gauss-radau', 'gauss-radau-right', 'equi','equi-noleft'};
min_nodes = [1, 1, 2, 2, 2, 2, 1];
a = 0.0;
b = 1.0;

for tt=1:length(type)
    
    for nn = min_nodes(tt):max_nodes(tt)
        
        % Create collocation object
        coll = collocation(a, b, nn, type{tt});
              
        % The 'evaluate_integral' routine of the collocation class operates
        % on a cell array of objects of the solution class, using the
        % provided axpy function. Here, it is checked that for a scalar
        % solution, the result computes to straightword evaluation of the
        % collocation formula as a scalar product weights x point values.
        y    = zeros(nn,1);
        ysol = cell(1,nn);
        for kk=1:nn
            y(kk,1)    = randn(1);
            ysol{1,kk} = solution_arraylinear(y(kk,1), 0.0, 1);
        end
        
        coll = coll.set_node_solutions(ysol);
        int  = coll.evaluate_integral;
        assert( abs(int.y - dot(coll.weights, y))<1e-12, 'Function coll.evaluate_integral failed to give same result as dot(weights, y)');
        
        % For every pp < order, the collocation formula should integrate
        % the polynomial f(x) = x^pp exactly.
        for pp=1:coll.order-1
            
            % Create cell array of scalar_solution objects, one for each
            % node x(kk) ; the value at x_i is x(kk)^pp
            fy = cell(1,nn);
            for kk=1:nn
                
                % Create a scalar_solution with value equal to kk-th nodes
                % and right hand side function f(x) = x^pp
                temp     = solution_arraynonlinear(coll.nodes(kk), @(x) x.^pp, @(x) pp*x.^(pp-1), 1.0);
                fy{1,kk} = temp.f;
            end
            
            % Set solution at node x_i to f(x_i)
            coll = coll.set_node_solutions(fy);
            
            % Now evaluate integral by collocation
            int = coll.evaluate_integral;
           
            % Anti-derivative of x^pp is (1/pp+1)*x^(pp+1).
            int_ex = (1/(pp+1))*(b^(pp+1) - a^(pp+1));
            
            % Make sure the polynomial x^pp is integrated exactly
            assert( abs(int.y - int_ex) < 1e-12, '%s : Failed to integrate polynomial of order %2i exactly with method of order %2i', type{tt}, pp, coll.order);
        end
    end

   % fprintf('Test successful : %s \n', type{tt});

end

% Test that the get_end_value routine computes u0 + sum q_j*f(u_j)
for tt=1:length(type)
   for M=2:10
      
       lambda = -0.66;
       mass   = 1.0;
       u0     = 1.0;
       sol_0  = solution_arraylinear(u0, lambda, mass);
       coll   = collocation(0.0, 1.0, M, type{tt});
       coll   = coll.set_initial_solution(sol_0);
       
       u      = randn(M,1);
       u_sols = cell(1,M);
       for mm=1:M
          u_sols{1,mm} = solution_arraylinear(u(mm,1), lambda, mass); 
       end
       coll = coll.set_node_solutions(u_sols);
       
       u_end         = u0 + dot(coll.weights, lambda*u);
       sol_u_end_mat = solution_arraylinear(u_end, lambda, mass);
       
       sol_u_end     = coll.get_end_value();
       diff          = sol_u_end.axpy(-1.0, sol_u_end_mat);
       assert( diff.getnorm < 1e-14, 'Function get_end_value does not return the expected result');
       
   end
end
fprintf('[0] - COLLOCATION - All tests successful \n');