classdef solution_arraynonlinear < solution
    
    % This class represents a solution of the nonlinear initial value problem
    %
    % y' = f(y)
    %
    % with y f(y) being linear arrays. Also requires a function handle Df,
    % that evalutes the Jacobian of f and is used in IMP_SOLVE for a Newton
    % method.
    
    properties(SetAccess=private,GetAccess=public)
        y;
        rhs;
        Df;
        M;
    end
    
    methods(Access=public)
        
        %% CONSTRUCTOR
        function obj = solution_arraynonlinear(y, f, Df, M)
            % Constructor function.
            
            % Sanity checks
            assert( size(y,2)==1, 'Solution has to be an array');
            assert( isa(f, 'function_handle'), 'f has to be a function handle');
            assert( isa(Df, 'function_handle'), 'Df has to be a function handle');
            % Make sure both f and Df can be evaluated on y and that Df(y)
            % is of size length(y) x length(y)
            assert( length(f(y))==size(Df(y),1), 'Df(y) must return a matrix with a number of rows equal to the length of f(y)');
            assert( length(f(y))==size(Df(y),2), 'Df(y) must return a matrix with a number of columns equal to the length of f(y)');
                        
            obj.y   = y;
            obj.rhs = f;
            obj.Df  = Df;
            obj.M   = M;
        end
        
        %% FUNCTIONS IMPLEMENTING THE ABSTRACT FUNCTIONS FROM SUPERCLASS
        
        function obj = axpy(obj, a, y)
            % Computes a*x + y
            
            assert( length(a)==1, 'Input a has to be a scalar');
            assert( isa(y, 'solution_arraynonlinear'), 'Input y for function addsol belonging to an object of type solution_arraynonlinear also must be of type solution_arraynonlinear.');
            
            obj.y = a*obj.y + y.y;
        end
        
        function obj = f(obj)
            % Evaluate f for current solution and update it to f(y).
            
            obj.y = obj.rhs(obj.y);
        end
        
        function obj = applyM(obj)
            obj.y = obj.M*obj.y;
        end
        
        function obj = imp_solve(obj,alpha,varargin)
            % Solve the nonlinear problem
            % y_sol - alpha*f(y_sol) = y
            % for y_sol using Newton's method with a fixed tolerance of
            % 1e-14 and a maximum of 100 iterations.
            
            assert( length(alpha)==1, 'Input alpha has to be a scalar');
            
            incr = ones(length(obj.y),1);
            x_n0 = obj.y;
            iter = 0;
            
            % Solve the nonlinear system 
            % y - alpha*f(y) = y_alt <=> y - alpha*f(y) - y_alt = 0
            % using Newton's method.
            while norm(incr) > 1e-14
               Df_n = obj.M - alpha*obj.Df(x_n0); % Jacobian
               f_n  = obj.M*x_n0 - alpha*obj.rhs(x_n0) - obj.y; % f
               % Solve linear system for increment
               incr = Df_n\(-f_n);
               x_n0 = x_n0 + incr;
               iter = iter+1;
               
               assert( iter<=250, 'Newton iteration did not converge in 250 iterations, increment inf last iteration = %5.3e', norm(incr));
            end
            
            %fprintf('Newton method terminated with norm of increment %5.2e \n', norm(incr));
            
            % Overwrite solution with solution of nonlinear system
            obj.y = x_n0;
            
        end
        
        function obj = imp_solve_stages(obj, A, c, dt, t0, y0)
           error('solution_arraynonlinear:NotImplemented','The function imp_solve_stages for the class solution_arraynonlinear is not yet implemented'); 
        end
        
        %% GET  FUNCTIONS
        
        function nm = getnorm(obj)
            % Returns the norm of the current solution.
            
            nm = norm(obj.y);
        end
    end
end