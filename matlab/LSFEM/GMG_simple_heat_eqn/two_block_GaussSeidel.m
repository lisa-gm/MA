% ******************************************************************* %
% ****************** BLOCK GAUSS SEIDEL SIGMA & U ******************* %
% ******************************************************************* %

% we have J = [ J_ss, J_su; J_us, J_uu]

function sol = two_block_GaussSeidel(J, f, sol_init, max_iter, eps)
block_size = size(J,1)/2;

sigma = sol_init(1:block_size);
u = sol_init(block_size+1:end);

J_ss = J(1:block_size, 1:block_size);
J_su = J(1:block_size, block_size+1:end);
J_us = J(block_size+1:end, 1:block_size);
J_uu = J(block_size+1:end, block_size+1:end);

f_s = f(1:block_size);
f_u = f(block_size+1:end);

for it=1:max_iter
    sigma = sigma + J_ss \ ( f_s - [J_ss, J_su]*[sigma; u]);
    u = u + J_uu \ (f_u - [J_us, J_uu]*[sigma; u]);
    
    if(norm(J*[sigma; u] - f) < eps)
        fprintf('2 block GS converged after %d iterations, res: %d\n', it, norm(J*[sigma; u] - f));
        break;
    end
end

sol = [sigma; u];

end


