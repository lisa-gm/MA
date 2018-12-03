classdef collocation_sdc < collocation
    % In addition to normal collocation, collocation_sdc also computes and
    % stores the matrices Qmat and Smat
    
    properties(SetAccess=private,GetAccess=public)
        % The following two matrices are SDC specific: Assume that the
        % nodes for the collocation formulae are (t_m), m=0,...M
        %
        % Then, the rows of Smat correspond to approximations of the
        % integrals from t_m to t_m+1, that is given a vector of values
        % (u_m), m=0,...,M at the nodes t_m, the scalar product of the m-th
        % row of S with u is an approximation of the integral over u from
        % t_m to t_m+1 :
        %
        % I_m^m+1 := dot(S(m,:), u) ~= int_(t_m)^(t_m+1) u(t) dt
        %
        % The matrix-vector product Smat*u gives a vector in which the m-th
        % entry corresponds to I_m^m+1.
        %
        % The matrix Qmat corresponds to approximations of the integral
        % from xleft to t_m :
        %
        % I_0^m+1 := dot(Q(m,:), u) ~= int_(t_0)^(t_m) u(t) dt
        %
        Qmat
        Smat
               
        % Contains f(node_solutions), that is the f function of the
        % solution objects applied to all solutions in node_solutions
        F_node_solutions
        
        % Array containing the Delta_m
        delta_m
        
        % The Q_delta matrix build from the delta_m values
        Qdelta_mat
        
    end
    
    methods(Access=public)
        
        function obj = collocation_sdc(xleft,xright,nr_nodes, type)
            
            obj = obj@collocation(xleft, xright, nr_nodes, type);                                                
            obj = obj.buildSQmat;
            
            % Compute distances between nodes
            obj.delta_m    = zeros(1,nr_nodes);
            obj.Qdelta_mat = zeros(nr_nodes, nr_nodes);
            
            obj.delta_m(1,1) = obj.nodes(1,1) - xleft;
            obj.Qdelta_mat(1:nr_nodes,1) = obj.delta_m(1,1);
            for nn=2:obj.nr_nodes
                obj.delta_m(1,nn) = obj.nodes(1,nn)-obj.nodes(1,nn-1);
                obj.Qdelta_mat(nn:nr_nodes,nn) = obj.delta_m(1,nn);
            end
            
        end
        
        function obj = initialize(obj, y0)
            
            obj = obj.set_initial_solution(y0);
            % Initialize all node solutions with y0
            for nn=1:obj.nr_nodes
                obj.node_solutions{1,nn} = y0;
            end
            obj = obj.applyF;

        end
        
        function obj = applyF(obj)
            
            for nn=1:obj.nr_nodes
                obj.F_node_solutions{1,nn} = obj.node_solutions{1,nn}.f;
            end
        end
        
        function int = applySmat(obj, m)
            % Computes the scalar product Smat(m,:) x node_solutions and
            % returns the resulting solution object
            
            % Initiate with zero solution
            int = obj.F_node_solutions{1,1}.axpy(-1.0, obj.F_node_solutions{1,1});
            assert( int.getnorm==0, 'Tried to create a zero solution, but failed');
            
            for nn=1:obj.nr_nodes
                int = obj.F_node_solutions{1,nn}.axpy( obj.Smat(m,nn), int);
            end

        end
        
        function int = applyQmat(obj, m)
            % Computes the scalar product Smat(m,:) x node_solutions and
            % returns the resulting solution object
            
            % Initiate with zero solution
            int = obj.F_node_solutions{1,1}.axpy(-1.0, obj.F_node_solutions{1,1});
            
            assert( int.getnorm==0, 'Tried to create a zero solution, but failed');
            
            for nn=1:obj.nr_nodes
                int = obj.F_node_solutions{1,nn}.axpy( obj.Qmat(m,nn), int);
            end
            
        end
        
        function obj = applySweep(obj, theta, gamma)
            % Performs one sweep over all nodes
                      
            % Create zero solution
            zero = obj.F_node_solutions{1,1}.axpy(-1.0, obj.F_node_solutions{1,1});
            
            % First, update F_node_solutions
            obj = obj.applyF;
            
            % Initial step using y0
            int = obj.applySmat(1);  
                    % int <- sum s_(1,j)*f(U^k_j)
                    
            % Multiply with gamma
            int = int.axpy(gamma, zero);
                    % int <- gamma * sum s_(1,j)*f(U^k_j)
            
            int = obj.F_node_solutions{1,1}.axpy(-theta*obj.delta_m(1,1), int);
                    % int <- -theta*dt_1*f(U^k_1) + gamma * sum s_(1,j)*f(U^k_j)
                    
%             int = int.axpy(1.0, obj.y0);
            int = obj.y0.axpy(gamma, int);
                    % int <- gamma*y0 - theta*dt_1*f(U^k_1) + gamma * sum s_(1,j)*f(U^k_j)

            int = obj.node_solutions{1,1}.axpy(1-gamma, int);
                    % int <- gamma*y0 + (1-gamma)*U^k_1 - theta*dt_1*f(U^k_1) + gamma* sum s_(1,j) f(U^k_j)
                
            % In next step, we need the difference U^k_2 - U^k_1, so store in in auxiliary buffer
            diff_k_old = obj.node_solutions{1,1}.axpy(-1.0, obj.node_solutions{1,2});

            % then overwrite obj.node_solutions{1,1} with new iterate
            obj.node_solutions{1,1} = int.imp_solve( theta*obj.delta_m(1,1) , obj.node_solutions{1,1});   
                    % obj.node_solutions{1,1} <- U^(k+1)_(1) = gamma*y0 + (1-gamma)*U^k_1 + theta*dt_1* [ f(U^(k+1)_1) - f(U^k_1) ] + gamma*sum s_(1,j)*f(U^k_j)
            
            
            % Now perform the node-to-node sweep
            for nn=2:obj.nr_nodes
                
                int = obj.applySmat(nn);
                    % int <- sum s_(nn,j) f(U^k_j)
                    
                int = int.axpy(gamma, zero);
                    % int <- gamma * sum s_(nn,j) f(U^k_j)
                    
                int = obj.F_node_solutions{1,nn}.axpy(-theta*obj.delta_m(1,nn), int);
                    % int <- -theta*dt_nn*f(U^k_nn) + gamma * sum s_(nn,j) f(U^k_j)
                                        
                int = int.axpy(1.0, obj.node_solutions{1,nn-1});
                    % int <- U^(k+1)_(nn-1) - theta*dt_nn*f(U^k_nn) + gamma * sum s_(nn,j) f(U^k_j)

                int = diff_k_old.axpy(1-gamma, int);
                     % int <- U^(k+1)_(nn-1) + (1-gamma)*( U^k_nn - U^k_(nn-1) ) - theta*dt_m*f(U^k_nn) + gamma * sum s_(nn,j) f(U^k_j)
                     
                % Store difference U^k_(nn+1) - U^k_nn for next iteration
                if nn<obj.nr_nodes
                    diff_k_old = obj.node_solutions{1,nn}.axpy(-1.0, obj.node_solutions{1,nn+1});
                     % diff_k_old <- U^k_(nn+1) - U^k_nn
                end
                
                obj.node_solutions{1,nn} = int.imp_solve(theta*obj.delta_m(1,nn), obj.node_solutions{1,nn});
                    % obj.node_solutions{1,nn} <- U^(k+1)_nn = U^(k+1)_(nn-1) + (1-gamma)*( U^k_nn - U^k_(nn-1) ) + theta*dt_nn*[ f(U^(k+1)_nn) - f(U^k_nn) ] + gamma * sum s_(nn,j)*f(U^k_j)
            end

        end
        
        function res = get_residual(obj)
           % Compute the residual 
          
           obj = obj.applyF;
           
           res_v = zeros(obj.nr_nodes,1);
           u_v   = zeros(obj.nr_nodes,1);
           
           for nn=1:obj.nr_nodes
               
               int = obj.applyQmat(nn);      % -> Q(nn,:)*F(U)
               int = obj.y0.axpy(1.0, int);  % -> y0 + Q(nn,:)*F(U)
               int = obj.node_solutions{1,nn}.axpy(-1.0, int);
                
               res_v(nn,1) = int.getnorm;
               u_v(nn,1)   = obj.node_solutions{1,nn}.getnorm;
           end
           
           res = norm(res_v, inf)/norm(u_v, inf);
           
        end
        
        function [Mit, M0] = getScalarSweepMatrix(obj, lambda, theta, gamma)
            % For the linear test equation y' = lambda*y, a sweep of SDC(theta,gamma)
            % can be easily written as a matrix so that u_np1 = Mit*u_n0 + M0*u0
            Mit = ( speye(obj.nr_nodes) - theta*obj.Qdelta_mat*lambda ) \ ( (1-gamma)*speye(obj.nr_nodes) + ( gamma*obj.Qmat - theta*obj.Qdelta_mat)*lambda );
            M0  = ( speye(obj.nr_nodes) - theta*obj.Qdelta_mat*lambda ) \ ( gamma*speye(obj.nr_nodes) );
        end
    end
    
    methods(Access=private)
        
        
        function obj = buildSQmat(obj)
            % Builds first the matrix Qmat and then computes Smat from it
            
            % Allocate QMat. Note that M+1 = obj.get_nr_nodes.
            obj.Qmat = zeros(obj.nr_nodes, obj.nr_nodes);
            
            % To compute weights for approximation of the integrals from
            % xleft to some arbitrary point a, we need to compute integrals
            % from xleft to a over the Lagrangian polynomials corresponding
            % to obj.nodes
            poly = cell(1,obj.nr_nodes);
            for ii=1:obj.nr_nodes
                coeff       = zeros(1,obj.nr_nodes);
                coeff(1,ii) = 1.0;
                poly{1,ii}  = interp_poly(obj.nodes, coeff);
            end
            
            % For every node t_m, m=0,...,M, compute the (M+1) weights
            % q_(m,j), j=0,...,M for the quadrature rule approximating the
            % integral int_(xleft)^(t_m)
            for mm=1:obj.nr_nodes
                
                % Create a new collocation formula on the interval
                % [xleft, obj.nodes(mm)].
                % Because the Lagrangian polynomials to be integrated over
                % this interval are of order obj.get_nr_nodes-1, set the
                % number of nodes to the minimun value that still
                % guarantees exact integration
                coll_xleft_tm = collocation(obj.xleft, obj.nodes(mm), ceil(obj.nr_nodes/2), 'gauss-legendre');
                
                % The weight q_(m,j) is given by the integral
                % int_(xleft)^(t_m) over the j-th Lagrangian polynomial
                for jj=1:obj.nr_nodes
                    
                    % Evaluate jj-th Lagrangian polynomial on nodes of
                    % coll_xleft_tm
                    pi              =  poly{1,jj}.evaluate(coll_xleft_tm.nodes);
                    
                    % And evaluate quadrature to get weight q_(mm,jj)
                    obj.Qmat(mm,jj) = coll_xleft_tm.evaluate_integral_on_array(pi);
                end
                
                clear coll_xleft_tm
                
            end
            
            clear poly
            
            obj.Smat = collocation_sdc.buildSmat(obj.Qmat);
        end
        
    end
    
    methods(Access=private,Static=true)
        
        function Smat = buildSmat(Qmat)
            
            Smat = zeros(size(Qmat));
            Smat(1,:) = Qmat(1,:);
            for ii=2:size(Qmat,1)
                Smat(ii,:) = Qmat(ii,:) - Qmat(ii-1,:);
            end
            
        end
    end
end