classdef solution_arraylinear < solution
    
    % This class represents a solution of the linear initial value problem
    %
    % M*y' = A*y
    %
    
    properties(SetAccess=protected,GetAccess=public)
        y;
        A;
        M;
    end
    
    methods(Access=public)
        
        %% CONSTRUCTOR
        function obj = solution_arraylinear(y, A, M)
            % Constructor function.
            
            % y = value of the solution to set
            assert( size(y,2)==1, 'Solution has to be an array');
            assert( length(y)==size(A,1) && length(y)==size(A,2), 'A has to be a matrix of size length(y) x length(y)');
            assert( length(y)==size(M,1) && length(y)==size(M,2), 'M has to be a matrix of size length(y) x length(y)');
            obj.y = y;
            obj.A = A;
            obj.M = M;
            
        end
        
        %% FUNCTIONS IMPLEMENTING THE ABSTRACT FUNCTIONS FROM SUPERCLASS
        function obj = axpy(obj, a, y)
            % Computes a*x + y
            
            assert( length(a)==1, 'Input a has to be a scalar');
            assert( isa(y, 'solution_arraylinear'), 'Input y for function addsol belonging to an object of type solution_arraylinear also must be of type solution_arraylinear.');
            
            obj.y = a*obj.y + y.y;
        end
        
        function obj = f(obj)
            % Evaluate f for current solution and update it to f(y).
            
            obj.y = obj.A*obj.y;
        end
        
        function obj = applyM(obj)
            obj.y = obj.M*obj.y;
        end
        
        function obj = imp_solve(obj,alpha,varargin)
            % Solve the system (M - alpha*A)*y_sol = y and overwrite y
            % with y_sol.
            
            assert( length(alpha)==1, 'Input alpha has to be a scalar');
            
            LL = obj.M - alpha*obj.A;
            obj.y = LL\obj.y;
        end
        
        function stages = imp_solve_stages(obj, A, c, dt, t0, y0)
                
            % Solves the system of equations
                assert( isa(y0, 'solution_arraylinear'), 'Parameter y0 has to be of type solution_arraylinear');
                assert( length(t0)==1, 't0 has to be a scalar');
                assert( length(dt)==1, 'dt has to be a scalar');
                assert( length(c)==size(A,1) && length(c)==size(A,2), 'A has to be of size sxs and c of length s');
                
                s   = length(c);
                LL  = kron(speye(s),obj.M) - dt*kron(A,obj.A);
                rhs = dt*kron(A*ones(s,1), obj.A*y0.y );
                z   = LL\rhs;
                
                % Now create and return a cell array of solution objects
                % corresponding to the s stages
                stages = cell(1,s);
                for ii=1:s
                    inds         = (ii-1)*length(obj.y)+1:ii*length(obj.y);
                    stages{1,ii} = solution_arraylinear(z(inds), obj.A, obj.M);
                end
                
        end
                
        %% GET  FUNCTIONS
        function nm = getnorm(obj)
            % Returns the norm of the current solution.
            
            nm = norm(obj.y);
        end
    end
end