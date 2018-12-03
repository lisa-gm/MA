%
% Tests the summed collocation class
%
clearvars -except test_all_dir test_all_ll test_all_classname

addpath('../src');

T0       = 0.0;
T1       = 1.77;
Nsteps_v = 1:5;
types    = {'gauss-lobatto','gauss-legendre','gauss-radau','gauss-chebyshev','equi','equi-noleft'};
M_vec    = 2:5;

for tt=1:length(types)
    for mm=1:length(M_vec)
        
        err = zeros(1,length(Nsteps_v));
        
        for nn = 1:length(Nsteps_v)
            
            Tmesh  = linspace(T0, T1, Nsteps_v(nn)+1);
            
            sc = summed_collocation(Tmesh, M_vec(mm), types{tt});
            
            % compute values of a piecewise polynomial of order p-1
            u  = @(x) sin(pi*x);
            df = @(x) 0*x; % not needed
           
            % create cell arrays with solution objects
            uf = cell(1,M_vec(mm));
            for ii=1:Nsteps_v(nn)
                for jj=1:M_vec(mm)
                    uf{1,jj} = solution_arraynonlinear( u(sc.collocations{ii}.nodes(jj)), u, df, 1.0);
                end
                sc = sc.set_node_solutions(uf, ii);
            end
            
            int = sc.evaluate_integral;
            
            intex = solution_arraynonlinear( -(1/pi)*(cos(pi*T1) - cos(pi*T0)), u, df, 1.0);
            
            diff = int.axpy(-1.0, intex);
            
            err(1,nn) = diff.getnorm/intex.getnorm;
        end
        
        convrate = log(err(1,2:end)./err(1,1:end-1))./log(Nsteps_v(1:end-1)./Nsteps_v(2:end));
        assert( abs(min(convrate)-sc.collocations{1}.order)<1e-1 || min(convrate)>=sc.collocations{1}.order, 'Failed to vertify convergence order');
    end
end

warning('Need to implement test of summed_collocation sweep...')

fprintf('[0] - SUMMED_COLLOCATION - All tests successful \n');