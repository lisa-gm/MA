%%%%% HEAT EQN 2D %%%%%%%


% adapt RHS according to boundary conditions
N=10;
h = 1/N;
row_pts = (N-1);
total_pts = (N-1)*(N-1);
% grid, with x and y indices
[X, Y] = meshgrid(1/N:h:1-1/N, 1-1/N:-h:1/N);

f = sin(X.*Y);

% how to get boundary subtracted? 
% we have North, South, West, East bdy
% find corresponding indices
% corner gets two values subracted ...?

Nbdy = zeros(1, row_pts); 
% Sbdy = zeros(1, row_pts); 
Ebdy = linspace(0,1,row_pts)'; % from top to bottom 
Wbdy = linspace(0,1,row_pts)'; % i.e. north to south

% if we know just move everything up by 1

%Nbdy = ones(1,row_pts);
Sbdy = ones(1, row_pts); 
%Ebdy = ones(row_pts,1); % from top to bottom 
%Wbdy = ones(row_pts,1); % i.e. north to south

% have to take bdy into account
f(end, :) = f(end,:) + 1/h^2*Sbdy;
f(1, :) = f(1,:) + 1/h^2* Nbdy;

f(:,1) = f(:, 1) + 1/h^2*Wbdy;
f(:,end) = f(:,end) + 1/h^2* Ebdy;

L_h = set_up_L_h_2D(N);

% now transform everything into vectors
% to solve linear system of eqns

f_vec = zeros(total_pts,1);
for i=1:row_pts
    for j=1:row_pts
        f_vec((i-1)*row_pts+j) = f(row_pts+1-i, j);
    end 
end


% now finally solve L_h*u = f
u_vec = L_h \ f_vec;

% turn u into the right matrix format again
u = zeros(row_pts, row_pts);

for i=1:row_pts
    for j=1:row_pts
        u(row_pts+1-i, j) = u_vec((i-1)*row_pts+j);
    end
end 


% can i somehow check my solution later, i.e. 
%u_int=@(x,y) interp2d(X,Y,u,x,y); % and then compute delta(u)

[X_comp, Y_comp] =  meshgrid(0:h:1, 1:-h:0);

% compute corner pts, mean from other two
cp_NW = (Nbdy(1)+Wbdy(1))/2;
cp_NE = (Nbdy(end)+Ebdy(1))/2;
cp_SW = (Sbdy(1)+Wbdy(end))/2;
cp_SE = (Sbdy(end)+Ebdy(end))/2;

% assemble whole matrix
u_comp = [cp_NW, Nbdy, cp_NE; horzcat(Wbdy, u, Ebdy); cp_SW, Sbdy, cp_SE];
surf(X_comp,Y_comp, u_comp);
