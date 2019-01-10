% ******************************************************************** %
% ********** CONSTRUCT LIST WITH ADAPTIVELY REFINED MESH ************* %
% ******************************************************************** %

% want to keep quotient (d*ht)/(hx)^2 roughly 1
% ie d*ht = hx^2,
% assume that we go from ht_fine = 0.5*ht_coarse, then
% --> hx_fine = sqrt(d*0.5*ht_coarse), but we want nested meshes

% start with coarse grid pts in space, decide time accordingly
function [Nx_elem_list, Nt_elem_list] = adapt_refine_mesh(hx, Nx_elem, T, diff_const, levels)
S = Nx_elem*hx;
Nx_pts = Nx_elem + 1;

Nt_pts = ceil(diff_const*T/S^2*Nx_elem^2) + 1;

% now refine for different levels
Nx_elem_list = zeros(1, levels);
Nt_elem_list = zeros(1, levels);

Nx_elem_list(1) = Nx_pts - 1;
Nt_elem_list(1) = Nt_pts - 1;

for l=2:levels
    Nx_elem_list(l) = 2*Nx_elem_list(l-1);
    Nt_elem_list(l) = 4*Nt_elem_list(l-1);
end

Nx_elem_list = fliplr(Nx_elem_list);
Nt_elem_list = fliplr(Nt_elem_list);

end