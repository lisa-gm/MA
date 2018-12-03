% ***************************************************************** %
% ****** fct that returns indices of inner pts of u & sigma ******* %
% ***************************************************************** %

function [inner_ind_sigma, inner_ind_u] = inner_ind(bdy_cond, Nx_elem, Nt_elem)
Nx_pts = Nx_elem + 1;
Nt_pts = Nt_elem + 1;
tot_pts = Nx_pts*Nt_pts;

if(strcmp(bdy_cond, 'Dirichlet'))
    inner_ind_u = [];
    for ind=2:Nt_pts
        inner_ind_u = [inner_ind_u, (ind-1)*Nx_pts+2:ind*Nx_pts-1];
    end
    inner_ind_sigma = 1:tot_pts;
end

if(strcmp(bdy_cond, 'Neumann'))
    inner_ind_sigma = [];
    for ind=1:Nt_pts
        inner_ind_sigma = [inner_ind_sigma, (ind-1)*Nx_pts+2:ind*Nx_pts-1];
    end
    inner_ind_u = Nx_pts+1:tot_pts;
end

end