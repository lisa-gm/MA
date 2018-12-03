classdef solution_xv_lennardjones < solution_xv
    
    
    properties(SetAccess=private,GetAccess=public)
        x
        v
        Mx
        Mv
        epsilon
        sigma
    end
    
    methods(Access=public)
        
        function obj = solution_xv_lennardjones(x, v, Mx, Mv, epsilon, sigma)
            
            assert( length(x) > 1, 'Need at least two entries in x');
            assert( length(x)==length(v), 'Input x and v must be arrays of the same length');
            assert( length(x)==size(Mx,1) && length(x)==size(Mx,2), 'Input Mx must be of dimension lenght(x) x length(x)');
            assert( length(x)==size(Mv,1) && length(x)==size(Mv,2), 'Input My must be of dimension lenght(x) x length(x)');
            assert( length(epsilon)==1, 'Input epsi must be a scalar');
            assert( length(sigma)==1,   'Input sigma must be a scalar');
            assert( size(x,2)==1, 'Input x has to be a column vector');
            assert( size(v,2)==1, 'Input v has to be a column vector');
            
            obj.x  = x;
            obj.v  = v;
            obj.Mx = Mx;
            obj.Mv = Mv;
            obj.epsilon = epsilon;
            obj.sigma   = sigma;
            
        end
        
        function obj = axpy(obj, a, y)
            
            assert(length(a)==1, 'Input a must be a scalar');
            assert(isa(y,'solution_xv_lennardjones'), 'Input y must of an object of class solution_xv_lennardjones');
            
            obj.x = a*obj.x + y.x;
            obj.v = a*obj.v + y.v;
        end
        
        function obj = f1(obj)
            % It is x' = v, so f1(x,v) = v
            obj.x = obj.v;
            obj.v = 0*obj.v;
        end
        
        function obj = f2(obj)
            % It is v' = f2(x)
            obj.v = solution_xv_lennardjones.lennardjones(obj.x, obj.epsilon, obj.sigma);
            obj.x = 0*obj.x;
        end
        
        function obj = f(obj)
            % Evaluate (x,v)' = (v, f(x))
            vv    = solution_xv_lennardjones.lennardjones(obj.x, obj.epsilon, obj.sigma);
            obj.x = obj.v;
            obj.v = vv;
            
        end
        
        function obj = applyM(obj)
            obj.x = obj.Mx*obj.x;
            obj.v = obj.Mv*obj.v;
        end
        
        function obj = imp_solve(obj, alpha)
            % Solves x - alpha*f(x) = y using Broyden's method
            assert( length(alpha)==1, 'Input alpha has to be a scalar');
            
            % If alpha = 0, invert only the mass matrix which requires a
            % simple linear solve
            if alpha==0
                xneu = obj.Mx\obj.x;
                vneu = obj.Mv\obj.v;
            else
                % IMPLEMENTATION BASED ON MATLAB'S FSOLVE FUNCTION
                %                 opts = optimset('TolFun',1e-14,'TolX',1e-14,'TolCon',1e-14,'MaxIter', 500);%, 'Display', 'Off');
                %
                %                 n = length(obj.x);
                %
                %                 % Current x and v define the right hand side
                %                 b = [obj.x ; obj.v];
                %
                %                 ff  = @(x) x - alpha*[ x(n+1:2*n) ; solution_xv_lennardjones.lennardjones(x(1:n), obj.epsilon, obj.sigma) ] - b;
                %                 xx  = fsolve(ff, b, opts);
                %                 assert( norm(ff(xx), inf) < 1e-10, 'Computed solution is not a root. Residual = %5.3e \n', norm(ff(xx), inf));
                %
                %                 % Overwrite x and v with solution of nonlinear system
                %                 obj.x = xx(1:n);
                %                 obj.v = xx(n+1:2*n);
                
                % BROYDEN'S METHOD (its orders of magnitude faster then
                % fsolve
                n     = length(obj.x);
                res   = 1.0;
                
                xalt  = obj.x;
                valt  = obj.v;
                
                % (x,v) - alpha*f(x,v) - y
                Fxalt = xalt - alpha*valt - obj.x;
                Fvalt = solution_xv_lennardjones.lennardjones(xalt, obj.epsilon, obj.sigma);
                Fvalt = valt - alpha*Fvalt - obj.v;
                
                Jac   = speye(2*n);
                iter  = 0;
                while (res > 1e-12)
                    
                    % Compute increment
                    Deltaxv = Jac\( -[Fxalt; Fvalt] );
                    xneu    = xalt + Deltaxv(1:n);
                    vneu    = valt + Deltaxv(n+1:2*n);
                    
                    nmxv = norm(Deltaxv);
                    %fprintf('Residual in Broydens method : %5.3e \n', res);
                    
                    Fxneu = xneu - alpha*vneu - obj.x;
                    Fvneu = solution_xv_lennardjones.lennardjones(xneu, obj.epsilon, obj.sigma);
                    Fvneu = vneu - alpha*Fvneu - obj.v;
                    
                    DeltaF = [Fxneu;Fvneu] - [Fxalt;Fvalt];
                    
                    res = norm( [Fxneu;Fvneu], inf);
                    
                    % Update estimated Jacobian
                    Jac = Jac + (1/nmxv^2)*(DeltaF - Jac*Deltaxv)*Deltaxv.';
                    
                    % Prepare next iteration
                    xalt  = xneu;
                    valt  = vneu;
                    Fxalt = Fxneu;
                    Fvalt = Fvneu;
                    iter  = iter+1;
                    assert( iter<=200, 'Broydens method failed to converge after 200 iterations');
                end
            end
            
            obj.x = xneu;
            obj.v = vneu;
        end
        
        function obj = imp_solve_stages(obj, A, c, dt, t0, y0)
            error('solution_xv_lennardjones:NotImplemented','The method imp_solve_stages is not implemented for the class solution_xv_lennardjones');
        end
        
        function n = getnorm(obj)
            n = sqrt( norm(obj.x).^2 + norm(obj.v).^2 );
        end
        
        function e = get_kinetic_energy(obj)
            e = sum(obj.v.^2)*0.5;
        end
        
        function e = get_potential_energy(obj)
            e = 0.0;
            for nn=1:length(obj.x)
                for mm=nn+1:length(obj.x)
                    r = norm(obj.x(nn) - obj.x(mm));
                    f = 4*obj.epsilon*( (obj.sigma/r).^12 - (obj.sigma/r).^6 );
                    e = e + f;
                end
            end
        end
        
        function e = get_energy(obj)
            e = obj.get_kinetic_energy + obj.get_potential_energy;
        end
    end
    
    methods(Access=protected,Static=true)
        
        function forces = lennardjones(x,e,s)
            
            forces = zeros(length(x),1);
            
            for nn=1:length(x)
                
                % Position of particle one
                r1 = x(nn);
                
                % Note: f_(n,n) = 0
                for mm=1:nn-1
                    
                    % Position of particle two
                    r2 = x(mm);
                    r  = r1 - r2;
                    rn = norm(r);
                    
                    % evaluate interaction
                    f  = ((48*e)/(rn^2))*( (s/rn)^12 - 0.5*(s/rn)^6)*r;
                    
                    % Add contribution f_(n,m)
                    forces(nn,1) = forces(nn,1) + f;
                    
                    % Add contribution f_(m,n) = -f_(n,m)
                    forces(mm,1) = forces(mm,1) - f;
                end
                
            end
            
        end
        
    end
    
end