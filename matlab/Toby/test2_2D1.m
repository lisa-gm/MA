%clear all;
% close all;
clc;

%fem 2d test for debugging the 4D - diffusion works

%==============================
% parameters
%==============================

ele_dim = [40,80];                                                  %elems per dim

dof_dim = ele_dim + 1;                                              %dofs per dim

xmin = [0,0];                                                      %grid dims (x space, y time)
xmax = [+1,+2];
xrng = xmax - xmin;

ele_h = xrng./ele_dim;

%==============================
% storage
%==============================

%for plots
X = zeros(dof_dim(1),dof_dim(2));                                   %coord/soln for plot 
Y = zeros(dof_dim(1),dof_dim(2));

U1 = zeros(dof_dim(1),dof_dim(2));                                  %solutions for plot
U2 = zeros(dof_dim(1),dof_dim(2));
U3 = zeros(dof_dim(1),dof_dim(2));
U4 = zeros(dof_dim(1),dof_dim(2));

%working memory
bf_eval = zeros(4,1);                                               %basis evaluations 
bf_grad = zeros(4,2);

%quadrature
qpt_np = 3;
qpt_pp = 0.5*([-sqrt(3/5),0,+sqrt(3/5)]+1);                         %3-point gauss on [0,1]
qpt_ww = 0.5*[5/9,8/9,5/9];

%fem matrices
mtx_n = prod(dof_dim);                                              %problem size

F = zeros(mtx_n,1);                                                 %rhs
M = zeros(mtx_n);                                                   %mass
K = zeros(mtx_n);                                                   %stiffness

%==============================
% functions
%==============================


%bc
fn_bc1    	= @(x,y) sin(pi*x)+1;                                             
fn_bc_mask1	= @(x,y) (abs(y) == 0); 
%fn_bc_mask1	= @(x,y) (abs(x)<0.2)&&(y<0.5);    

%bc
fn_bc2    	= @(x,y) 1;                                             
fn_bc_mask2	= @(x,y) (abs(x)==xmax(1)) || (abs(x) == 0) ;   
 

%rhs
fn_f      	= @(x,y) pi^2*sin(pi*x)*cos(pi*y);                 	
fn_f_mask   = @(x,y) 1;                                                         

fn_pos2idx = @(i,j) (j-1)*dof_dim(1) + (i-1) + 1;                   %pos to idx 

fn_bf_eval = @(x,y) [   (1-x)*(1-y);                                %basis functions
                        (  x)*(1-y);
                        (1-x)*(  y);
                        (  x)*(  y)];
                
fn_bf_grad = @(x,y) [   (-1)*(1-y)/ele_h(1),(1-x)*(-1)/ele_h(2);
                        (+1)*(1-y)/ele_h(1),(  x)*(-1)/ele_h(2);
                        (-1)*(  y)/ele_h(1),(1-x)*(+1)/ele_h(2);
                        (+1)*(  y)/ele_h(1),(  x)*(+1)/ele_h(2)];

fn_dof_idx = @(i,j) [   fn_pos2idx(i,j);
                        fn_pos2idx(i+1,j);
                        fn_pos2idx(i,j+1);
                        fn_pos2idx(i+1,j+1)];

%==============================
% loop dofs
%==============================

counter = 0;

for dof_j = 1:dof_dim(2)

    for dof_i = 1:dof_dim(1)
        
%         fprintf('dof %2d %2d %2d\n',fn_pos2idx(dof_i,dof_j),dof_i,dof_j);
        
        X(dof_i,dof_j) = xmin(1) + (dof_i-1)*ele_h(1);                                                          %store coords
        Y(dof_i,dof_j) = xmin(2) + (dof_j-1)*ele_h(2);
        
        U1(dof_i,dof_j) = fn_bc_mask1(X(dof_i,dof_j),Y(dof_i,dof_j))*fn_bc1(X(dof_i,dof_j),Y(dof_i,dof_j));
        U2(dof_i,dof_j) = fn_bc_mask2(X(dof_i,dof_j),Y(dof_i,dof_j))*fn_bc2(X(dof_i,dof_j),Y(dof_i,dof_j));
%         U2(dof_i,dof_j) = fn_f_mask(X(dof_i,dof_j),Y(dof_i,dof_j))*fn_f(X(dof_i,dof_j),Y(dof_i,dof_j));         %calc data and bc
        

    end
    
end

%==============================
% loop elements
%==============================

%loop elements - L2 projection
for ele_j = 1:ele_dim(2)
    
    for ele_i = 1:ele_dim(1)
   
        ele_idx = fn_pos2idx(ele_i,ele_j);                                  %index for assembly
        
        ele_ref_x = xmin(1) + (ele_i-1)*ele_h(1);                           %reference vertex
        ele_ref_y = xmin(2) + (ele_j-1)*ele_h(2);
        
        dof_idx = fn_dof_idx(ele_i,ele_j);                                  %dof indices for assembly
        
        fprintf('ele	%3d %2d %2d | %f %f \n',ele_idx,ele_i,ele_j,ele_ref_x,ele_ref_y);
                
        %==============================
        % loop quad points
        %==============================
        
        for qpt_j = 1:qpt_np
        
            for qpt_i = 1:qpt_np                                        	%tensor qpts
   
                qpt_x = ele_ref_x + qpt_pp(qpt_i)*ele_h(1);                 %global coords
                qpt_y = ele_ref_y + qpt_pp(qpt_j)*ele_h(2);
                
                f_eval = fn_f_mask(qpt_x,qpt_y)*fn_f(qpt_x,qpt_y);       	%evaluate the data function 
                
                bf_eval = fn_bf_eval(qpt_pp(qpt_i),qpt_pp(qpt_j));          %evaluate basis (ref coords)
                bf_grad = fn_bf_grad(qpt_pp(qpt_i),qpt_pp(qpt_j));
                
                %fprintf('qpt	%3d %2d | %f %f | %f %f %+f\n',qpt_i,qpt_j,qpt_pp(qpt_i),qpt_pp(qpt_j),qpt_x,qpt_y,f_eval);

%                 %==============================
%                 % loop vtxs - evaluate basis
%                 %==============================
%                 
%                 for vtx_j = 0:1
% 
%                     for vtx_i = 0:1
% 
%                         idx = vtx_j*2+vtx_i+1;
%                       
%                         bf_eval(idx) = (1-vtx_i+(2*vtx_i-1)*qpt_pp(qpt_i))*(1-vtx_j+(2*vtx_j-1)*qpt_pp(qpt_j)); %eval on ref
%                         
%                         bf_grad(idx,:) = [(2*vtx_i-1)*(1-vtx_j+(2*vtx_j-1)*qpt_pp(qpt_j)),(1-vtx_i+(2*vtx_i-1).*qpt_pp(qpt_i))*(2*vtx_j-1)];
%                        
%                     end
% 
%                 end
                
                %==============================
                % tensor dofs - assembly
                %==============================
                
                for dof_j = 1:4    %loop bfs/rows/test

                    test_integral(dof_idx(dof_j)) =  test_integral(dof_idx(dof_j)) + 
                    F(dof_idx(dof_j)) = F(dof_idx(dof_j)) + bf_eval(dof_j)*f_eval*(qpt_ww(qpt_i)*qpt_ww(qpt_j));            %ass rhs                                                  %rhs
                    
                    for dof_i = 1:4                                                                                         %loop bfs/cols/trial
                        
                        
                        M(dof_idx(dof_i),dof_idx(dof_j)) = M(dof_idx(dof_i),dof_idx(dof_j)) + bf_eval(dof_i)*bf_eval(dof_j) * (qpt_ww(qpt_i)*qpt_ww(qpt_j)); %mass
                        
                        %K(dof_idx(dof_i),dof_idx(dof_j)) = K(dof_idx(dof_i),dof_idx(dof_j)) + bf_grad(dof_i,:)*bf_grad(dof_j,:)' * (qpt_ww(qpt_i)*qpt_ww(qpt_j));   %stiffness - diffusion
    
                        %spacetime also works see loop3
                        K(dof_idx(dof_i),dof_idx(dof_j)) = K(dof_idx(dof_i),dof_idx(dof_j)) + bf_grad(dof_i,1)*bf_grad(dof_j,1) * (qpt_ww(qpt_i)*qpt_ww(qpt_j));   %stiffness - 1x1 spacetime du/dx
                        K(dof_idx(dof_i),dof_idx(dof_j)) = K(dof_idx(dof_i),dof_idx(dof_j)) + bf_eval(dof_i)*bf_grad(dof_j,2) * (qpt_ww(qpt_i)*qpt_ww(qpt_j));   %stiffness - 1x1 spacetime du/dt
                        
                    end 

                end
                
            end
    
        end
        
    end
    
end

%==============================
% boundary condition
%==============================

%loop dofs
for dof_j = 1:dof_dim(2)

    for dof_i = 1:dof_dim(1)
        
        if fn_bc_mask1(X(dof_i,dof_j),Y(dof_i,dof_j))                %test mask
   
            row_idx = fn_pos2idx(dof_i,dof_j);                    	%select dof row

            K(row_idx,:) = zeros(1,prod(dof_dim));                  %force K to identity
            K(row_idx,row_idx) = 1;

            F(row_idx) = fn_bc1(X(dof_i,dof_j),Y(dof_i,dof_j));   	%set soln to bc value
        
        end
        
        if fn_bc_mask2(X(dof_i,dof_j),Y(dof_i,dof_j))                %test mask
   
            row_idx = fn_pos2idx(dof_i,dof_j);                    	%select dof row

            K(row_idx,:) = zeros(1,prod(dof_dim));                  %force K to identity
            K(row_idx,row_idx) = 1;

            F(row_idx) = fn_bc2(X(dof_i,dof_j),Y(dof_i,dof_j));   	%set soln to bc value
        
        end   
    
    end

end

save('diff_opt.mat', 'K');
  
%==============================
% solve
%==============================

U = K\F;

%F = M*F;        %show rhs in basis

%==============================
% loop dofs
%==============================


for dof_j = 1:dof_dim(2)
   
    for dof_i = 1:dof_dim(1)
        
        U3(dof_i,dof_j) = F(fn_pos2idx(dof_i,dof_j));
        U4(dof_i,dof_j) = U(fn_pos2idx(dof_i,dof_j));               %write solution to mesh
      
    end
    
end

%==============================
% disp
%==============================

% subplot(1,4,1);
% surf(X,Y,U1);
% 
% subplot(1,4,2);
% surf(X,Y,U2);
% 
% subplot(1,4,3);
% surf(X,Y,U3);
% 
% subplot(1,4,4);
% surf(X,Y,U4);

figure
surf(X,Y,U4);
title('Approx Solution');
xlabel('space');
ylabel('time');

%cond(K)    
    