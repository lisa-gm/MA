classdef solution
    
    % The abstract class "Solution" defines the general properties that are 
    % required from and object representing a solution y of the initial 
    % value problem
    %
    % M*y' = f(y)
    %
    % A subclass representing a specific type of solution has to implement 
    % the following functionality
    % 
    % 1. "addsol"    : Add two solutions 
    % 2. "multsol"   : Multiply a solution with a scalar
    % 3. "f"         : Evaluate f and return f(y)
    % 4. "imp_solve" : Solve the system (Id - alpha*f)(y_sol) = y
    %                  and return y_sol
    % 5. "getnorm"   : Return a norm ||y|| of y
    
    methods(Abstract=true)
                
        % Computes a*obj + y where a is a scalar and y is another object of
        % the solution class (typically, y also has to be of the same
        % subclass as obj).
        axpy(obj, a, y)
        
        % Evaluate the right hand side f for solution y and overwrite y
        % with f(y)
        % NOTE: Here, for the time being, the assumption is that y and f(y)
        % are objects of the same type, e.g. vectors in R^n. This need not
        % be so and more complex applications will likely require a more
        % general approach here.
        f(obj)
        
        % Overwrites the solution y with M*y
        applyM(obj)
        
        % Solve the linear or nonlinear (depending on f) problem
        % (M - alpha*f)(y_sol) = y
        % and overwrite y with the solution y_sol.
        % NOTE: INTERFACE PROBABLY HAS TO BE EXTENDED TO ACCOMODATE
        % TIME-DEPENDENT RIGHT HAND SIDES AND MASS MATRICES!
        imp_solve(obj, alpha, varargin)
   
        % In order to use a solution object within a higher order implicit
        % Runge-Kutta method, where a large system for s stages has to be
        % solved, it has to also provide the following routine which
        % returns the solution of a system of s equations
        %
        % M*z_i = dt*sum(j=1..s) a(i,j)*f(y0 + z_j, t0 + c_j*dt)
        %
        % for i=1,...,s.
        imp_solve_stages(obj, A, c, dt, t0, y0)
        
        % Return a norm of the solution.
        getnorm(obj)
                
    end
end