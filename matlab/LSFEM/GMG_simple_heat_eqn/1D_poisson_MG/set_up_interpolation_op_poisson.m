% ******************************************************************** %
% *************** CONSTRUCT INTERPOLATION OPERATOR ******************* %
% ******************************************************************** %

% ONLY CHANGED INTERPOLATION OPERATOR

% Nx_elem_list goes from coarse to fine
function [I, R] = set_up_interpolation_op_poisson(Nx_elem_list, levels)
I = {};
R = {};

%  first build INTERPOLATION OPERATOR for either u or sigma

for l=1:(levels-1)
Nx_pts_c = Nx_elem_list(l+1)+1;
Nx_pts = Nx_elem_list(l)+1;

temp = ones(Nx_pts,1);
mat_temp = spdiags([temp, 2*temp, temp], -1:1, Nx_pts, Nx_pts);
%R{j} = 1/2*mat_temp(2:2:end, :);
I_one = (1/2*mat_temp(1:2:end, :))';

% need this operator twice, for sigma and u
Int = zeros(2*size(I_one));
Int(1:Nx_pts, 1:Nx_pts_c) = I_one;
Int(Nx_pts+1:end, Nx_pts_c+1:end) = I_one;

if(l == 1)
    bdy_ind = get_bdy_ind_poisson(Nx_pts);
    Int(bdy_ind, :) = 0;
end

I{l} = Int;
R{l} = Int';
end

end