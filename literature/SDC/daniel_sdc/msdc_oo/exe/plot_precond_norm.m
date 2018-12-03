clear all;
close all;
beep off;

addpath('../src');
addpath('../src/quad_nodes');
addpath('../../mparareal_oo/src');

Tend   = 1.0;
lambda = -10 + 0i;
nodes_v = 3:1:15;

norm_qdelta = zeros(1,length(nodes_v));
norm_esti   = zeros(1,length(nodes_v));

for nn=1:length(nodes_v)
    nodes = nodes_v(nn);
    coll    = collocation_sdc(0, Tend, nodes, 'gauss-lobatto');
    
    Qdelta = coll.Qdelta_mat;
    Sdelta = diag(diag(Qdelta));
    LLmat = tril(ones(nodes));
    
    assert( norm(Qdelta - LLmat*Sdelta, inf)<1e-14, 'Failed: Q_d = L*S_d');
    
    LLmat_inv = spdiags([-ones(nodes,1) ones(nodes,1)], -1:0, nodes, nodes);
    LLmat_inv = full(LLmat_inv);
    
    assert( norm(speye(nodes) - LLmat_inv*LLmat, inf) < 1e-14, 'Failed: Inverse of L');
    
    cc    = 1 - lambda*coll.delta_m;
    assert( min(abs(cc))>= 1, 'Found c_i with absolute value smaller than one');
    
    CCmat = spdiags([-ones(nodes,1) cc.'], -1:0, nodes, nodes);
    CCmat = full(CCmat);
    
    % build inverse of CCmat from formula
    CCmat_inv = zeros(nodes);
    for jj=1:nodes
        CCmat_inv(jj,jj) = 1./cc(jj);
        for ii=jj+1:nodes
            CCmat_inv(ii,jj) = (1./cc(ii))*CCmat_inv(ii-1,jj);
        end
    end
    assert( norm(inv(CCmat) - CCmat_inv, inf) < 1e-14, 'Failed: Formula for inverse of C');
    
    Q_precond = speye(nodes) - lambda*coll.Qdelta_mat;
    
    % build inverse of Id - lambda*Q_delta = L*C
    Q_precond_inv = zeros(nodes);
    for jj=1:nodes
        Q_precond_inv(jj,jj) = 1/cc(jj);
        
        if jj<nodes
            Q_precond_inv(jj+1,jj) = (1/cc(jj+1))*(1/cc(jj) - 1);
            
            for ii=(jj+2):nodes
                Q_precond_inv(ii,jj) = (1/cc(ii))*Q_precond_inv(ii-1,jj);
            end
            
        end
        
    end
    assert( norm( inv(Q_precond) - Q_precond_inv, inf)< 1e-14, 'Failed: Formulate for inverse of preconditioner');
       
    % Verify that (Id - lambda*Qdelta) = L*( inv(L) - lambda*Sdelta )
    assert( norm(Q_precond - LLmat*CCmat, inf) < 1e-14, 'Failed: Id - lambda*Qdelta = L*CC');
    
    norm_qdelta(1,nn) = norm(inv(Q_precond), inf);
    temp = 0;
    for mm=1:nodes
        temp = temp+prod(1./cc(mm:nodes));
    end
    norm_esti(1,nn) = sum( abs(1./abs(cc(1:end-1)) - 1.0 )  ) + 1/abs(cc(end));
    %norm_esti(1,nn) = 1./(1 - abs(lambda));
    
    clf;
    plot(nodes_v, norm_esti, 'r-'); hold on;
    plot(nodes_v, norm_qdelta, 'ro--', 'markerfacecolor', 'r');
end