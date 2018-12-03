clear all;
close all;
beep off;

lw = 1.0;
fs = 12;
ms = 10;

addpath('../src');
addpath('../src/quad_nodes');

% Solves y' = lambda*y with SDC(theta,gamma) iteration
lambda  = -10 + 0*1i;
Tend    = 1.0;
uex     = exp(lambda*Tend);

Nsteps_v = 1:5;
nodes    = 4;

type = 'gauss-lobatto';
theta = 1.0;
Nit   = 7;

gamma = ones(1,Nit);
% gamma(1) = 0.5;
% gamma(2) = 0.75;

err   = zeros(1,length(Nsteps_v));

for ll=1:length(Nsteps_v)
    
    Tmesh = linspace(0, Tend, Nsteps_v(ll)+1);
    u_ini = 1.0;
    sc = summed_collocation(Tmesh, nodes, type);
    
    for jj=1:Nsteps_v(ll)
                
        un0 = u_ini*ones(nodes,1);
        for ii=1:Nit
            [Mit, M0] = sc.collocations{jj}.getScalarSweepMatrix(lambda, theta, gamma(ii));
            unp = Mit*un0 + M0*(u_ini*ones(nodes,1));
            un0 = unp;
        end
        
        % prepare next timestep
       	u_ini = un0(nodes,1);
    end
    
    err(1,ll) = abs(un0(nodes,1) - uex);
end

fprintf('\n');
for nn=1:length(Nsteps_v)
   fprintf('Error: %5.2e \n', err(1,nn)); 
end

convrate = log10(err(1,2:end)./err(1,1:end-1))./log10(Nsteps_v(1,1:end-1)./Nsteps_v(1,2:end));
fprintf('\n');
for nn=1:length(Nsteps_v)-1
   fprintf('Convergence rate: %5.2f \n', convrate(nn)); 
end