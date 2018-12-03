classdef solution_arraylinear_iter < solution_arraylinear
    
    properties(SetAccess=private,GetAccess=public)
        xleft;
        xright;
        N;
        nu;
        method;
        Nit;
    end
    
    methods(Access=public)
        
        %% CONSTRUCTOR
        function obj = solution_arraylinear_iter(y, A, M, xleft, xright, nu, method, Nit)
            
            % Constructor function.
            obj = obj@solution_arraylinear(y, A, M);
            
            % y = value of the solution to set
            obj.xleft = xleft;
            obj.xright = xright;
            obj.N = size(A,1);
            obj.nu = nu;
            obj.method = method;
            obj.Nit = Nit;
            
        end
        
        % Overload the imp_solve function
        function obj = imp_solve(obj,alpha,varargin)
            
            [M, M_b, bla, blub] = getTwoGridMat(obj.xleft, obj.xright, obj.N, obj.method, 1.0, -alpha*obj.nu);
                        
            assert( isa(varargin{1},'solution'), 'Initial guess has to be a solution object');
            
            rhs = obj.y;
            un0 = varargin{1}.y;
            %un0 = obj.y;
            for it=1:obj.Nit
                un1 = M*un0 + M_b*rhs;
                un0 = un1;
            end
            obj.y = un1;
           % res = norm( rhs - (un1 - alpha*obj.A*un1), inf);
           % fprintf('Final residual %5.3e \n', res);
            
        end
    end
end