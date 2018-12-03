classdef solution_heateq < solution
    
    % Describes a solution of the initial value problem
    %
    % y' = f(y)
    %
    % where y is a vector of point values and f corresponds to a finite
    % difference discretization of the 1D Laplacian with periodic boundary
    % conditions.
    
    properties(GetAccess=public,SetAccess=private)
        y     % solution
        Dxx   % discrete Laplacian
        Id    % identity
        xmesh % spatial mesh
    end
    
    methods(Access=public)
        
        %% CONSTRUCTOR
        function obj = solution_heateq(xleft, xright, y)
            % Constructor function.
            
            % The number of spatial mesh points is determined by the length
            % of the initial value y.
            Nx = length(y);
            assert( (xleft <= xright), 'Left endpoint xl has to be smaller than right endpoint xr');
            
            xaxis = linspace(xleft,xright,Nx+1);
            
            % remove last boundary: Because of periodic BC, its identical to first
            % point
            obj.xmesh = xaxis(1:Nx);
            
            % Now build discrete Laplacian
            dx      = abs(xaxis(2)-xaxis(1));
            e       = ones(Nx+1,1);
            obj.Dxx = spdiags([e -2*e e], -1:1, Nx, Nx);
            
            % Incorporate periodic BC
            obj.Dxx(1,Nx) = 1;
            obj.Dxx(Nx,1) = 1;
            
            % Scale with mesh width
            obj.Dxx = (1/dx^2)*obj.Dxx;
            
            % Create sparse identity matrix
            obj.Id  = speye(size(obj.Dxx));
            
            obj.y   = y;
            
        end
        
        %% FUNCTIONS IMPLEMENTING THE ABSTRACT FUNCTIONS FROM SUPERCLASS
        
        function obj = axpy(obj, a, y)
            % Computes a*x + y
            
            assert( length(a)==1, 'Input a has to be a scalar');
            assert( isa(y, 'solution_heateq'), 'Input y for function addsol belonging to an object of type solution_heateq also must be of type solution_heateq.');
            
            obj.y = a*obj.y + y.y;
        end
        
        function obj = f(obj)
            % Evaluates f by multiyplying y with Dxx and overwrites y with
            % f(y)
            
            obj.y = obj.Dxx*obj.y;
        end
        
        function obj = applyM(obj)
           % do nothing 
        end
        
        function obj = imp_solve(obj,alpha)
            % Solves the linear system (Id - alpha*Dxx)*y_sol = y and
            % overwrites y with the solution y_sol
            
            obj.y = (obj.Id - alpha*obj.Dxx)\obj.y;
        end
        
        function obj = imp_solve_stages(obj, A, c, dt, t0, y0)
           error('solution_heateq:NotImplemented','The function imp_solve_stages for the class solution_heateq is not yet implemented'); 
        end
        
        %% GET FUNCTIONS
        
        function y = getsol(obj)
            % Return value of solution
            y = obj.y;
        end
        
        function xmesh = getmesh(obj)
            % Return underlying spatial mesh
            xmesh = obj.xmesh;
        end
        
        function dx = getdx(obj)
            % Return meshsize of underlying spatial mesh
            dx = obj.xmesh(2) - obj.xmesh(1);
        end
        
        function Nx = getnx(obj)
            % Return number of degrees of freedom
            Nx = length(obj.xmesh);
        end
        
        function nm = getnorm(obj)
            % Return norm of solution
            nm = norm(obj.y, inf);
        end
    end
end