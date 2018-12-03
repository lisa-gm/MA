classdef collocation
    
    properties(GetAccess=public,SetAccess=protected)
        
        % The type of the collocation nodes
        type
        
        % Collocation nodes
        nodes
        
        % Number of collocation nodes
        nr_nodes
        
        % The collocation weights
        weights
        
        % Left endpoint of interval
        xleft
        
        % Right endpoint of interval
        xright
        
        % Order of collocation formulat, depends on type and nr_nodes
        order
        
        % Lebesgue constant for given nodes
        lebesgue
        
    end
    
    properties(GetAccess=public,SetAccess=protected)
        
        % A prescribed solution at xleft
        y0
        
        % A list of solutions, one at each node
        node_solutions
        
    end
    
    methods(Access=public)
        
        function obj = collocation(xleft,xright,nr_nodes,type)
            
            assert( xleft <= xright, 'Left endpoint must be smaller then right endpoint.');
            obj.xleft  = xleft;
            obj.xright = xright;
            obj.nr_nodes = nr_nodes;
            
            switch type
                
                case {'gauss-legendre', 'g-leg','Gauss-Legendre'}
                    assert( nr_nodes >=1, 'For Gauss-Legendre nodes, need nr_nodes >= 1');
                    
                    obj.type     = type;
                    [obj.nodes, obj.weights] = GaussLegendre(nr_nodes, xleft, xright);
                    obj.order    = 2*obj.nr_nodes;
                    
                case {'gauss-lobatto','g-lob','Gauss-Lobatto'}
                    assert( nr_nodes >=2, 'For Gauss-Lobatto nodes, need nr_nodes >= 2');
                    
                    obj.type     = type;
                    [obj.nodes, obj.weights] = GaussLobatto(nr_nodes, xleft, xright);
                    obj.order    = 2*obj.nr_nodes-2;
                    
                case {'gauss-radau','g-rad','Gauss-Radau'}
                    assert( nr_nodes >= 2, 'For Gauss-Radau nodes, need nr_nodes >= 2');
                    obj.type = type;
                    [obj.nodes, obj.weights] = GaussRadau(nr_nodes, xleft, xright);
                    
                    obj.order = 2*obj.nr_nodes-1;
                    
                case {'gauss-radau-right', 'g-rad-right', 'Gauss-Radau-Right'}
                    assert( nr_nodes >= 2, 'For Gauss-Radau nodes, need nr_nodes >= 2');
                    obj.type    = type;
                    obj.nodes   = GaussRadau_RightPoint(nr_nodes, xleft, xright);   
                    obj.weights = collocation.getWeights(obj.nodes, xleft, xright);
                    obj.order   = 2*obj.nr_nodes-1;
                    
                case {'gauss-chebyshev', 'g-cheb','Gauss-Chebyshev'}
                    obj.type = type;
                    [obj.nodes, obj.weights] = obj.Chebyshev(nr_nodes, xleft, xright);
                    % Because so far this is not used for weighted
                    % integrals, the order is just n
                    obj.order = obj.nr_nodes;
                    
                case {'gauss-perturbed'}
                    
                    obj.type     = type;
                    [obj.nodes, dummy] = GaussLegendre(nr_nodes, xleft, xright);
                    obj.nodes(1:5) = 1.1*obj.nodes(1:5);
                    obj.nodes(6:end) = 0.8*obj.nodes(6:end);
                    obj.weights = collocation.getWeights(obj.nodes, xleft, xright);
                    obj.order    = 2*obj.nr_nodes;
                    
                case {'equi'}
                    
                    obj.type = type;
                    obj.nodes = (linspace(xleft, xright, nr_nodes));
                    obj.weights = collocation.getWeights(obj.nodes, xleft, xright);
                    obj.order = obj.nr_nodes;
                    
                case {'equi-noleft'}
                    obj.type = type;
                    nodes_temp = linspace(xleft, xright, nr_nodes+1);
                    obj.nodes = nodes_temp(2:nr_nodes+1);
                    obj.weights = collocation.getWeights(obj.nodes, xleft, xright);
                    obj.order = obj.nr_nodes;
                    
                otherwise
                    error('collocation:NodeType','No implementation for selected node type');
            end
            
            % Allocate node_solutions and F_node_solutions
            obj.node_solutions   = cell(1,obj.nr_nodes);
            
            obj.lebesgue = lebesgue_constant(nr_nodes, obj.nodes, linspace(xleft, xright, 1e4));
            
            assert(size(obj.weights,1)==1, 'Weights has to be a row vector');
            assert(size(obj.nodes,1)==1, 'Nodes has to be a row vector');
            if obj.nr_nodes>=2
                % If nodes coincide, difference can be negative but only of the order of machine precision
                assert( min(diff(obj.nodes))>=-1e-12, 'Node vector has to be monotonically increasing');
            end
        end
        
        function obj = set_initial_solution(obj, y0)
            assert( isa(y0, 'solution'), 'Parameter y0 has to be an object of class solution');
            obj.y0 = y0;
        end
        
        function obj = set_node_solutions(obj, node_solutions)
            %
            
            assert( iscell(node_solutions), 'node_solutions has to be a cell array of solution class objects');
            assert( size(node_solutions,1)==1 && size(node_solutions, 2)==obj.nr_nodes, 'node_solutions must be a cell row array with nr_nodes many entries');
            
            for nn=1:obj.nr_nodes
                assert( isa(node_solutions{1,nn}, 'solution') , 'All entries in node_solutions must be objects of type solution');
                obj.node_solutions{1,nn} = node_solutions{1,nn};
            end
        end
        
        function int = evaluate_integral(obj)
            % Computes the sum over weights(i)*node_solutions(i).
            
            % Initialize with zero, using axpy with x=y and a=-1
            int = obj.node_solutions{1,1}.axpy(-1.0, obj.node_solutions{1,1});
            assert( int.getnorm==0, 'Tried to create a zero element, but failed');
            
            for nn=1:obj.nr_nodes;
                int = obj.node_solutions{1,nn}.axpy(obj.weights(nn), int);
            end
            
        end
        
        function int = get_end_value(obj)
            % Updates y_end
            
            % Initialize with zero, using axpy with x=y and a=-1
            int = obj.node_solutions{1,1}.axpy(-1.0, obj.node_solutions{1,1});
            assert( int.getnorm==0, 'Tried to create a zero element, but failed');
            
            for nn=1:obj.nr_nodes
                temp = obj.node_solutions{1,nn}.f;
                int  = temp.axpy(obj.weights(nn), int);
            end
            
            int = int.axpy(1.0, obj.y0);
            
        end
        
    end
    
    methods(Access=protected, Static=true)
        
        function [x,w] = Chebyshev(N, a, b)
            % A function to compute nodes and weights for Gauss-Chebyshev
            % quadrature.
            % --> This normally is for weighted integrals with a weighting
            % function w(x) = 1/sqrt(1-x^2). For non-weighted integrals, it
            % will only be of order N.
            
            % First, get the corresponding nodes
            x = ChebyshevRoots( N, 'Tn', [a b]);
            
            % Then, compute the corresponding weights as integrals over
            % Lagragian polynomials
            w = collocation.getWeights(x, a, b);
            
            % reshape into column vector
            x = x;
            w = w;
        end
        
    end
    
    methods(Access=protected)
        
        function int = evaluate_integral_on_array(obj, fx)
            % This method evaluates the collocation formulate for an array
            % of numbers, instead on a cell array of solution class
            % objects. It is used only by the collocation class itself, to
            % compute weights for some arbitrary set of nodes.
            assert( length(fx)==obj.nr_nodes, 'Length of fx must be equal to number of nodes');
            int = dot(obj.weights, fx);
            
        end
    end
    
    methods(Access=private, Static=true)
        
        function w = getWeights(x, a, b)
            % For a general set of collocation nodes, the corresponding weights
            % can be retrieved by computing the integrals int_a^b over the
            % corresponding Lagrange polynomials. This is not very efficient,
            % though.
            
            assert( min(x) >= a, 'Left endpoint a has to be equal or smaller than all nodes in x');
            assert( max(x) <= b, 'Right endpoint b has to be equal or larger then all nodes in x');
            
            w = zeros(1,length(x));
            
            coll = collocation(a, b, length(x), 'Gauss-Legendre');
            
            for nn=1:length(x)
                coeff       = zeros(length(x),1);
                coeff(nn,1) = 1.0;
                poly        = interp_poly(x, coeff);
                
                % Evaluate the Lagrangian polynomial at the nodes of the
                % collocation formula represented by coll
                y           = poly.evaluate(coll.nodes);
                
                w(1,nn)     = coll.evaluate_integral_on_array(y);
            end
        end
        
    end
end