% ******************************************************************** %
% *************** CONSTRUCT INTERPOLATION OPERATOR ******************* %
% ******************************************************************** %


% Nx_elem_list goes from coarse to fine
function [I, R] = set_up_interpolation_op_strategic(Nx_elem_list, Nt_elem_list, bdy_cond, levels)
I = {};
R = {};

%  first build INTERPOLATION OPERATOR for either u or sigma

for l=1:(levels-1)
Nx_pts_c = Nx_elem_list(l+1)+1;
Nt_pts_c = Nt_elem_list(l+1)+1;

Nx_pts = Nx_elem_list(l)+1;
Nt_pts = Nt_elem_list(l)+1;

if(Nt_pts_c == 1)
    fprintf('INVALID GRID SIZE');
    return;
end

pos_2id_f = @(i,j) (Nx_elem_list(l)+1)*(j-1) + i;
pos_2id_c = @(i,j) (Nx_elem_list(l+1)+1)*(j-1) + i;

tot_pts_c = Nx_pts_c*Nt_pts_c;
tot_pts = Nx_pts*Nt_pts;

I_one = zeros(tot_pts, tot_pts_c);


for j=1:Nt_pts
   
   for i=1:Nx_pts
            % odd points in space so: 1,3,5, ..
       if(mod(i,2)==1)

                % coarse grid pts in space   
                if(mod(j-1,Nt_pts_c-1)==0)
                    ind_j = (j-1)/(Nt_pts_c - 1) + 1;
                    ind_i = (i+1)/2;
                    
                    I_one(pos_2id_f(i,j),pos_2id_c(ind_i, ind_j))=1;

                % fine grid pts in between, distribute weights dep on dist
                % to the two closest coarse grid pts
                else
                   ind_j = floor((j-1)/(Nt_pts_c - 1)) + 1;
                   ind_i = (i+1)/2;
                   
                   b = mod(j-1, Nt_pts_c-1)/(Nt_pts_c - 1);
                   a = 1 - b;

                   I_one(pos_2id_f(i,j),pos_2id_c(ind_i, ind_j))=a;
                   I_one(pos_2id_f(i,j),pos_2id_c(ind_i, ind_j+1))=b;
                end
             
            else        % meaning even points in space, ie 2,4,6, ...
              fprintf('i even\n'); 
            if(mod(j-1, Nt_pts_c - 1) == 0)
               
                fprintf('i : %d, j : %d\n', i, j);
                ind_j = (j-1)/(Nt_pts_c - 1) + 1;
                
                % 0.5 from left and right
                I_one(pos_2id_f(i,j),pos_2id_c(i/2, ind_j))=0.5;
                I_one(pos_2id_f(i,j),pos_2id_c(i/2+1, ind_j))=0.5;
                pos_2id_c(i/2, ind_j)
                pos_2id_c(i/2+1, ind_j)
            else
                fprintf('otherwise\n');
                ind_j = floor((j-1)/(Nt_pts_c - 1)) + 1;
                ind_j
                i/2
               pos_2id_c(i/2,ind_j)
                % from below and above, and left and right
                b = mod(j-1,Nt_pts_c-1)/(Nt_pts_c - 1);
                a = (1-b);
                a = a/2;
                b = b/2;
                
                % where to write them
                I_one(pos_2id_f(i,j),pos_2id_c(i/2,ind_j))=a;
                I_one(pos_2id_f(i,j),pos_2id_c(i/2+1,ind_j))=a;
                I_one(pos_2id_f(i,j),pos_2id_c(i/2,ind_j+1))=b;
                I_one(pos_2id_f(i,j),pos_2id_c(i/2+1,ind_j+1))=b;
            
            end              
       end    
    end
end

% need this operator twice, for sigma and u
Int = zeros(2*size(I_one));
Int(1:tot_pts, 1:tot_pts_c) = I_one;
Int(tot_pts+1:end, tot_pts_c+1:end) = I_one;

if(l == 1)
    bdy_ind = get_bdy_ind(Nx_pts, Nt_pts, bdy_cond);
    Int(bdy_ind, :) = 0;
end

I{l} = Int;
R{l} = Int';
end

end
