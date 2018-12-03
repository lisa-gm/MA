classdef interp_poly
    % An interpolation polynomial
    
    properties(GetAccess=public,SetAccess=private)
        
        % The nodes at which the interpolation conditions are specified
        nodes
        
        % The function values prescribed at the nodes
        ynodes
        
        % The coefficients of the resulting interpolation polynomial in the
        % Newton basis
        coeff
    end
    
    methods(Access=public)
        
        function obj = interp_poly(x,y)
            % Creates an instance of the class representing an interpolation polynomial through (x(i),y(i))
            
            assert( length(x)==length(y), 'Vectors x and y have to have the same length');
            obj.nodes  = x;
            obj.ynodes = y;
            
            % Compute coefficients of interpolation polynomial in Newton
            % representation
            obj.coeff = interp_poly.newpoly(x,y); 
            
        end
        
        function fyi = evaluate(obj, xi)
            
            % Horner scheme
            fyi = obj.coeff(obj.get_nr_nodes);
            for ii=2:obj.get_nr_nodes
               fyi = obj.coeff(obj.get_nr_nodes-ii+1) + (xi - obj.nodes(obj.get_nr_nodes-ii+1)).*fyi;
            end
        end
        
        function nr_nodes = get_nr_nodes(obj)
           nr_nodes = length(obj.nodes); 
        end
    end
    
    methods(Access=private,Static=true)
    
        function D = newpoly(X,Y)
            %---------------------------------------------------------------------------
            %NEWPOLY   Construction of the collocation polynomial.
            %          The method is Newton interpolation.
            % Sample calls
            %   [C] = newpoly(X,Y)
            %   [C,D] = lnewpoly(X,Y)
            % Inputs
            %   X   vector of abscissas
            %   Y   vector of ordinates
            % Return
            %   C   coefficient list for the Newton polynomial
            %   D   divided difference table
            %
            % NUMERICAL METHODS: MATLAB Programs, (c) John H. Mathews 1995
            % To accompany the text:
            % NUMERICAL METHODS for Mathematics, Science and Engineering, 2nd Ed, 1992
            % Prentice Hall, Englewood Cliffs, New Jersey, 07632, U.S.A.
            % Prentice Hall, Inc.; USA, Canada, Mexico ISBN 0-13-624990-6
            % Prentice Hall, International Editions:   ISBN 0-13-625047-5
            % This free software is compliments of the author.
            % E-mail address:      in%"mathews@fullerton.edu"
            %
            % Algorithm 4.5 (Newton Interpolation Polynomial).
            % Section	4.4, Newton Polynomials, Page 234
            %---------------------------------------------------------------------------
            
            n = length(X);
            D = zeros(n,n);
            D(:,1) = Y';
            for j=2:n,
                for k=j:n,
                    D(k,j) = (D(k,j-1)-D(k-1,j-1))/(X(k)-X(k-j+1));
                end
            end
            
            % Do not need coefficients for monomial basis here...
            
%             C = D(n,n);
%             for k=(n-1):-1:1,
%                 C = conv(C,poly(X(k)));
%                 m = length(C);
%                 C(m) = C(m) + D(k,k);
%             end

            % Required coefficients are on the diagonal of divided
            % differences tableau
            D = diag(D);
        end
        
    end
end