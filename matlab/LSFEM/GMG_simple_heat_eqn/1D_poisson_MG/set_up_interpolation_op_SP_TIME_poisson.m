% ******************************************************************** %
% *************** CONSTRUCT INTERPOLATION OPERATOR ******************* %
% ******************************************************************** %

% ONLY CHANGED INTERPOLATION OPERATOR

% Nx_elem_list goes from coarse to fine
function [I, R] = set_up_interpolation_op_SP_TIME_poisson(Nx_elem_list, Nt_elem_list, bdy_cond, levels)
I = {};
R = {};

%  first build INTERPOLATION OPERATOR for either u or sigma

for l=1:(levels-1)
Nx_pts_c = Nx_elem_list(l+1)+1;
Nt_pts_c = Nt_elem_list(l+1)+1;

Nx_pts = Nx_elem_list(l)+1;
Nt_pts = Nt_elem_list(l)+1;

tot_pts_c = Nx_pts_c*Nt_pts_c;
tot_pts = Nx_pts*Nt_pts;

I_one = zeros(tot_pts, tot_pts_c);

for j=1:Nt_pts
   
    % 1st, 3rd, 5th, ... row in time
    if(mod(j,2)==1)
        for i=1:Nx_pts

            % 1st, 3rd, 5th, ... point in space   
            if(mod(i,2)==1)
                I_one(i+(j-1)*Nx_pts,(i-1)/2+1+(j-1)/2*Nx_pts_c)=1;

            % 2nd, 4th, 6th,... point in space 
            else
            I_one(i+(j-1)*Nx_pts,[i/2+(j-1)/2*Nx_pts_c,i/2+1+(j-1)/2*Nx_pts_c])=[0.5,0.5];
            end
            
        end
       
    % 2nd, 4th, 6th, ... row in space
    else
        for i=1:Nx_pts

            % 1st, 3rd, 5th, ... point in space   
            if(mod(i,2)==1)
            I_one(i+(j-1)*Nx_pts,(i-1)/2+1+(j-2)/2*Nx_pts_c)=0.5;
            I_one(i+(j-1)*Nx_pts,(i-1)/2+1+(j)/2*Nx_pts_c)=0.5;

            % 2nd, 4th, 6th,... point in space 
            else
            I_one(i+(j-1)*Nx_pts,i/2+(j-2)/2*Nx_pts_c)=0.25;
            I_one(i+(j-1)*Nx_pts,i/2+(j)/2*Nx_pts_c)=0.25;
            I_one(i+(j-1)*Nx_pts,i/2+1+(j-2)/2*Nx_pts_c)=0.25;
            I_one(i+(j-1)*Nx_pts,i/2+1+(j)/2*Nx_pts_c)=0.25; 
            end           
        end    
    end
end

% need this operator twice, for sigma and u
Int = zeros(2*size(I_one));
Int(1:tot_pts, 1:tot_pts_c) = I_one;
Int(tot_pts+1:end, tot_pts_c+1:end) = I_one;

if(l == 1)
    bdy_ind = get_bdy_ind_poisson(Nx_pts, Nt_pts, bdy_cond);
    Int(bdy_ind, :) = 0;
end

I{l} = Int;
R{l} = Int';
end

end