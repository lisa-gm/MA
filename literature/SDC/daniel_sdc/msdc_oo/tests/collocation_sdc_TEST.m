%
% Tests the collocation_sdc class
%
clearvars -except test_all_dir test_all_ll test_all_classname

addpath('../src');
addpath('../src/quad_nodes');

% define function and its anti-derivative
f = @(x) erf(5*x);
T0      = -0.15;
T1      = 0.1;

types   = {'gauss-legendre','gauss-lobatto', 'gauss-radau', 'gauss-chebyshev', 'gauss-radau-right'};

% Test that the matrices Q and S satisfy the expect summation property:
% Q(m,:) = S(1,:) + ... + S(m,:)
for ii=1:length(types)
    for M=2:10
        y0 = solution_arraylinear(0.0, 1.0, 1.0);
        coll = collocation_sdc(T0, T1, M, types{ii});
        coll = coll.initialize(y0);
        Q    = coll.Qmat;
        S    = coll.Smat;
        
        for mm=1:M
            Ssum = sum(S(1:mm,:), 1);
            assert( norm(Q(mm,:) - Ssum)<1e-12, 'Q and S do not satisfy the expected summation property');
        end
    end
end

for ii = 1:length(types)
    
    % For too small node numbers, the collocation formulas can produce very
    % wrong results
    for nodes=[5 7 10 15]
        
        y0   = solution_arraynonlinear(0.0, f, @(x) 0*x, 1.0);
        coll = collocation_sdc(T0, T1, nodes, types{ii});
        coll = coll.initialize(y0);
        
        % Generate a cell array of solution objects, each representing the
        % solution at one collocation node
        node_solutions = cell(1,nodes);
        for nn=1:nodes
            node_solutions{1,nn} = solution_arraynonlinear(coll.nodes(nn), f, @(x) 0*x, 1.0);
        end
        % Put cell array of node solutions into coll
        coll = coll.set_node_solutions(node_solutions);
        % And update F_node_solutions
        coll = coll.applyF;
        
        % First test: Make sure that the result from the partial
        % collocation formulas give a somewhat reasonable approximation of
        % actual integrals (here computed with Matlab's quad function)
        for mm=1:nodes
            int_m      = coll.applyQmat(mm);
            int_incr_m = coll.applySmat(mm);
            
            int_m_ex      = quad(f, T0, coll.nodes(mm), 1e-12);
            if (mm==1)
                int_incr_m_ex = int_m_ex;
            else
                int_incr_m_ex = quad(f, coll.nodes(mm-1), coll.nodes(mm), 1e-12);
            end
            
            err_inc = abs(int_incr_m_ex - int_incr_m.y);
            err     = abs(int_m_ex      - int_m.y);
            % The tolerance of 1e-4 is quite arbitrary.... but it makes
            % sure that the actions of Qmat and Smat are somewhat close to
            % the integrals they are supposed to approximate
            assert( err_inc < 1e-4, '%s, nodes = %2i, mm = %2i : Error %5.2e in Smat probably too large', types{ii}, nodes, mm, err_inc);
            assert( err < 1e-4, '%s, nodes = %2i, mm = %2i : Error %5.2e in Qmat probably too large', types{ii}, nodes, mm, err);
        end
        
        % Second test: All the partial formulas should have an order equal
        % to the number of nodes, because the weights are computed as
        % integrals over Lagrangian polynomals. Make sure that polynomials
        % of degree up to nodes-1 are integrated exactly.
        for pp=1:nodes-1
            
            % Generate solution objects for every node
            for nn=1:nodes
                node_solutions{1,nn} = solution_arraynonlinear(coll.nodes(nn), @(x) x.^pp, @(x) pp*x.^(pp-1), 1.0);
            end
            coll = coll.set_node_solutions(node_solutions);
            coll = coll.applyF;
            
            for mm=1:nodes
                int_m         = coll.applyQmat(mm);
                int_incr_m    = coll.applySmat(mm);
                
                int_m_ex      = (1/(pp+1))*(coll.nodes(mm)^(pp+1) - T0^(pp+1));
                if (mm==1)
                    int_incr_m_ex = int_m_ex;
                else
                    int_incr_m_ex = (1/(pp+1))*(coll.nodes(mm)^(pp+1)-coll.nodes(mm-1)^(pp+1));
                end
                
                err_inc = abs(int_incr_m_ex - int_incr_m.y);
                err     = abs(int_m_ex      - int_m.y);
                assert( err < 1e-12, '%s, nodes = %2i, mm = %2i : Failed to integrate polynomial of order %2i exactly (in Qmat)', types{ii}, nodes, mm, pp);
                assert( err_inc < 1e-12, '%s, nodes = %2i, mm = %2i : Failed to integrate polynomial of order %2i exactly (in Smat)', types{ii}, nodes, mm, pp);
            end
        end
    end
    
    %fprintf('%s successful \n', types{ii});
end

%types   = {'gauss-legendre','gauss-lobatto', 'gauss-radau', 'gauss-chebyshev'};

% Test that the sweep in collocation_sdc corresponds to U^(k+1) = U_0 +
% Q*lambda*U^k in matrix form
for tt=1:length(types)
    for M=2:10
        
        lambda    = -0.66;
        mass      = 1.0;
        u_np1     = zeros(M,1);
        
        %
        % FIRST: Make sure that SDC with theta=0 corresponds to a matrix-wise Picard iteration
        %
        
        % node values
        u = randn(M,1);
        
        if strcmp(types{tt},'gauss-lobatto')
            % first node values is also initial value
            u0 = u(1);
        else
            % random initial value
            u0 = randn(1);
        end
        sol_0 = solution_arraylinear(u0, lambda, mass);
        
        % create collocation object
        coll = collocation_sdc(0.0, 0.77, M, types{tt});
        coll = coll.initialize(sol_0);
        
        % create a cell array of M solution objects
        node_sols = cell(1,M);
        for mm=1:M
            node_sols{1,mm} = solution_arraylinear(u(mm,1), lambda, mass);
        end
        
        % FIRST: Test Picard sweep
        
        % Set them in coll
        coll = coll.set_node_solutions(node_sols);
        
        % Now update them performing a sweep with theta=0, gamma=1 (Picard)
        coll = coll.applySweep(0.0, 1.0);
        
        % and rebuild a vector of nodes values
        for mm=1:M
            u_np1(mm,1) = coll.node_solutions{1,mm}.y;
        end
        
        % Now directly compute u_np1 using the Q matrix
        [Mit, M0] = coll.getScalarSweepMatrix(lambda, 0.0, 1.0);
        u_np1_mat = Mit*u + M0*(u0*ones(M,1));
        
        % and check that results are identical
        assert( norm(u_np1 - u_np1_mat) < 1e-12, 'Component-wise Picard sweep is not identical to matrix Picard iteration');
        
        %
        % SECOND: Make sure that an SDC sweep with theta=1 and gamma=1 is identical to its matrix formulation
        %
        
        % Reset node solutions
        coll = coll.set_node_solutions(node_sols);
        
        % Perform a sweep with theta=1 (SDC)
        coll = coll.applySweep(1.0, 1.0);
        
        % build node value vector
        for mm=1:M
            u_np1(mm,1) = coll.node_solutions{1,mm}.y;
        end
        
        % Compute directly
        [Mit, M0] = coll.getScalarSweepMatrix(lambda, 1.0, 1.0);
        u_np1_mat = Mit*u + M0*(u0*ones(M,1));
        
        assert( norm(u_np1 - u_np1_mat) < 1e-12, 'Component-wise SDC sweep is not identical to matrix SDC sweep');
        
        % THIRD: Test SDC(theta) sweep for random theta
        
        % perform 5 tests with random numbers
        for tests=1:5
            theta = randn(1);
            
            % Reset node solutions
            coll = coll.set_node_solutions(node_sols);
            
            % Perform a sweep with random value of theta
            coll = coll.applySweep(theta, 1.0);
            
            % build node value vector
            for mm=1:M
                u_np1(mm,1) = coll.node_solutions{1,mm}.y;
            end
            
            % Compute directly
            [Mit, M0] = coll.getScalarSweepMatrix(lambda, theta, 1.0);
            u_np1_mat = Mit*u + M0*(u0*ones(M,1));
            
            assert( norm(u_np1 - u_np1_mat) < 1e-12, 'Component-wise SDC(thetea) sweep is not identical to matrix SDC(theta) sweep');
        end
        
        %
        % THIRD: Test three random values for theta and make sure that applySweep is equivalant to the matrix formulation
        %
        
        % perform 5 tests, each with different random numbers
        for tests=1:5
            
            theta = randn(1);
            gamma = randn(1);
            
            % Reset node solutions
            coll = coll.set_node_solutions(node_sols);
            
            % Perform a sweep with random values for theta and gamma
            coll = coll.applySweep(theta, gamma);
            
            % build node vector
            for mm=1:M
                u_np1(mm,1) = coll.node_solutions{1,mm}.y;
            end
            
            % get iteration matrix for SDC(theta,gamma)
            [Mit, M0] = coll.getScalarSweepMatrix(lambda, theta, gamma);
            u_np1_mat = Mit*u + M0*(u0*ones(M,1));
            
            assert( norm(u_np1 - u_np1_mat) < 1e-12, 'Component-wise SDC(thete,gamma) sweep is not identical to matrix SDC(theta,gamma) sweep');
            
        end
    end
end


fprintf('[0] - COLLOCATION_SDC - All tests successful \n');