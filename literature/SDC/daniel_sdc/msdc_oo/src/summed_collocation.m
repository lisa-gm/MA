classdef summed_collocation
    
    properties(SetAccess=protected,GetAccess=public)
        
        % A list of collocation class objects
        collocations
        
        % The length of the individual time steps
        meshwidth
        
        % Number of total quadrature nodes in all time steps
        nr_nodes_total
        
        % array collecting all nodes
        nodes_total
        
        % List of indices: For every collocation object there is a list of
        % integers corresponding to the position of the nodes of this time
        % step if all nodes are put into one single array
        node_indices
        
        % Multi time step quadrature weight matrix
        Qmat
    end
    
    methods(Access=public)
        
        function obj = summed_collocation(Tmesh, nr_nodes, type)
            
            obj.nr_nodes_total = 0;
            obj.nodes_total    = zeros(1,(length(Tmesh)-1)*nr_nodes);
            
            for ii=1:length(Tmesh)-1
                obj.collocations{ii} = collocation_sdc(Tmesh(ii), Tmesh(ii+1), nr_nodes, type);
                obj.node_indices{ii} = obj.nr_nodes_total+1:obj.nr_nodes_total+nr_nodes;
                obj.nr_nodes_total   = obj.nr_nodes_total + nr_nodes;
                obj.nodes_total(obj.node_indices{ii}) = obj.collocations{ii}.nodes;
            end
            
            % Build diagonal of Qmat
            obj.Qmat = blkdiag(obj.collocations{1}.Qmat);
            for ii=2:length(Tmesh)-1
                obj.Qmat = blkdiag(obj.Qmat, obj.collocations{ii}.Qmat);
            end
            
            % For every time step except the last, build the full Q quadrature weight
            % matrix
            Qmats_full = cell(length(Tmesh)-2,1);
            for jj=1:length(Qmats_full)
                Qmats_full{jj,1} = zeros(obj.collocations{jj}.nr_nodes);
                for ii=1:obj.collocations{jj}.nr_nodes
                    Qmat_full{jj,1}(ii,:) = obj.collocations{jj}.weights.';
                end
            end
            
            for bb=1:length(Qmats_full)
                % fill all blocks below diagonal with Qmat_full{bb,1}
                for ii=bb+1:(length(Tmesh)-1)
                    obj.Qmat(obj.node_indices{ii},obj.node_indices{bb}) = Qmat_full{bb,1};
                end
            end
        end
        
        function obj = set_node_solutions(obj, node_solutions, index)
            
            assert( ( (1<=index) && (index<=length(obj.collocations))), 'Index must be between 1 and number of timesteps');
            assert( iscell(node_solutions), 'node_solutions has to be a cell array of solution class objects');
            assert( size(node_solutions,1)==1 && size(node_solutions, 2)==obj.collocations{index}.nr_nodes, 'node_solutions must be a cell row array with nr_nodes many entries');
            
            obj.collocations{index} = obj.collocations{index}.set_node_solutions(node_solutions);
        end
        
        function obj = initialize(obj, y0, index)
            assert( ( (1<=index) && (index<=length(obj.collocations))), 'Index must be between 1 and number of timesteps');
            obj.collocations{index} = obj.collocations{index}.initialize(y0);
        end
        
        function obj = applySweep(obj, theta, gamma, index)
            assert( ( (1<=index) && (index<=length(obj.collocations))), 'Index must be between 1 and number of timesteps');
            obj.collocations{index} = obj.collocations{index}.applySweep(theta, gamma);
        end
        
        function obj = applyF(obj)
           
            for nn=1:obj.get_nr_collocations
               obj.collocations{nn} = obj.collocations{nn}.applyF;
            end
        end
        
        function int = evaluate_integral(obj)
            
            % Initialize with zero, using axpy with x=y and a=-1
            int = obj.collocations{1}.node_solutions{1,1}.axpy(-1.0, obj.collocations{1}.node_solutions{1,1});
            assert( int.getnorm==0, 'Tried to create a zero element, but failed');
            
            for ii=1:obj.get_nr_collocations
                int_loc = obj.collocations{ii}.evaluate_integral;
                int = int_loc.axpy(1.0, int);
            end
            
        end
        
        function nr_colls = get_nr_collocations(obj)
            nr_colls = length(obj.collocations);
        end
    end
    
    methods(Access=private)
        % blablub
    end
end