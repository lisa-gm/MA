function bdy_ind = get_bdy_ind(Nx_pts, Nt_pts, bdy_cond)
tot_pts = Nx_pts*Nt_pts;
%for each column that corresponds to a bdy entry

if(strcmp(bdy_cond, 'Dirichlet'))
    bdy_ind_sigma = [];
    bdy_ind_u = 1:Nx_pts;
    for ts=2:Nt_pts
        bdy_ind_u = [bdy_ind_u,(ts-1)*Nx_pts+1, ts*Nx_pts];
    end
end

if(strcmp(bdy_cond, 'Neumann'))
    bdy_ind_u = 1:Nx_pts;
    bdy_ind_sigma = [];
    for ts=1:Nt_pts
        bdy_ind_sigma = [bdy_ind_sigma, (ts-1)*Nx_pts+1, ts*Nx_pts];
    end
end

bdy_ind = [bdy_ind_sigma, tot_pts + bdy_ind_u];

end
