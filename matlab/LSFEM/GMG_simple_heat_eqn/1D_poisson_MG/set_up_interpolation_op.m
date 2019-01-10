% ******************************************************************** %
% *************** CONSTRUCT INTERPOLATION OPERATOR ******************* %
% ******************************************************************** %


% Nx_elem_list goes from coarse to fine

% ONLY COARSEN IN SPACE (take out every other point in space, but leave 
% points in time) 

% stencil for inner values of the form
% 
%      1/8   -----   1/8
%      0.25 -- x -- 0.25
%      1/8   -----  1/8 

% boundary in space: 0.5 -- x -- 0.5 

function [I, R, Nx_elem_list, Nt_elem_list] = set_up_interpolation_op_SP(Nx_elem_list, Nt_elem_list, bdy_cond, levels)
I = {};
R = {};

%  first build INTERPOLATION OPERATOR for either u or sigma

for l=1:(levels-1)
Nx_pts_c = Nx_elem_list(l+1)+1;
Nt_pts_c = Nt_elem_list(l)+1;

Nx_pts = Nx_elem_list(l)+1;
Nt_pts = Nt_elem_list(l)+1;

% since we don't coarsen in time!
Nt_elem_list(l+1) = Nt_elem_list(l);

tot_pts_c = Nx_pts_c*Nt_pts_c;
tot_pts = Nx_pts*Nt_pts;

I_one = zeros(tot_pts, tot_pts_c);

% first put in values for the boundary 
for i=1:Nx_pts
    if(mod(i,2)==0)
        I_one(i, i/2) = 0.5;
        I_one(i, i/2+1) = 0.5;

        I_one(tot_pts-Nx_pts+i, tot_pts_c-Nx_pts_c+i/2) = 0.5;
        I_one(tot_pts-Nx_pts+i, tot_pts_c-Nx_pts_c+i/2+1) = 0.5;

    else
        I_one(i, (i-1)/2+1) = 1;
        I_one(tot_pts-Nx_pts+i, tot_pts_c-Nx_pts_c+(i-1)/2+1) = 1;

    end
end

for j=2:Nt_pts-1
        for i=1:Nx_pts
            
            if(mod(i,2)==1)
            % coarse grid points
            I_one(i+(j-1)*Nx_pts,(i-1)/2+1+(j-1)/2*Nx_pts_c)=1;

            % fine grid points
            else
            I_one(i+(j-1)*Nx_pts,[i/2+(j-1)/2*Nx_pts_c,i/2+1+(j-1)/2*Nx_pts_c])=[0.25,0.25];
            I_one(i+(j-1)*Nx_pts,i/2+(j-2)/2*Nx_pts_c)=0.125;
            I_one(i+(j-1)*Nx_pts,i/2+(j)/2*Nx_pts_c)=0.125;
            I_one(i+(j-1)*Nx_pts,i/2+1+(j-2)/2*Nx_pts_c)=0.125;
            I_one(i+(j-1)*Nx_pts,i/2+1+(j)/2*Nx_pts_c)=0.125; 
            end           
        end
end

% only need it once here, not LSFEM
Int = I_one;

if(l == 1)
    bdy_ind = get_bdy_ind(Nx_pts, Nt_pts, bdy_cond);
    Int(bdy_ind, :) = 0;
end

I{l} = Int;
R{l} = Int';
end

end