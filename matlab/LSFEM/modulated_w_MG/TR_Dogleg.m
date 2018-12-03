% ************************************************************** %
% ************* ADAPTED TRUST REGION DOGLEG ******************** %
% ************************************************************** %

function p = TR_Dogleg(grad_f, hess_f, update, TR_rad)
% determine tau

k = grad_f'*hess_f*grad_f;

%pb = - hess_f \ grad_f;
pb = update;

    if norm(pb) <= TR_rad
        % we can solve quad prob exactly 
        p = pb;
        
       %fprintf('newton step\n');
        
    else
        pu = -grad_f'*grad_f/k*grad_f;
        np_u = norm(pu);
        if np_u >= TR_rad
            
            % Use steepest descent direction
            p = TR_rad*pu/np_u;
            %fprintf('steepest descent\n');

        else
            % Find tau
            pc = pb - pu;
            coeffs = [norm(pc)^2, 2*pc'*pu, np_u^2-TR_rad^2];
            tau = max(roots(coeffs));
            p = pu + tau*pc;
            %fprintf('combined step\n');

        end
    end
end