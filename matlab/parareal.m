%%%%%% easy example from 50 years of time parallel time
% page 12 %%%%% lions, maday, turinci

% consider y' = -ay on [0,T], y(0)=y0, i.e. f(t,y) = a*y



%%%%%%%%%%% coarse approximation %%%%%%%%%%%%%%%
a = -1;
T = 5;
tsc = 5;           % time steps coarse
t_h = T/tsc;        % step size

max_iter = 3;
delta0 = zeros(tsc+1,1);
eps = 0.01;

error = zeros(1,max_iter+1);

% set up time vector t
t = linspace(0,T, tsc+1);

% initial value y_0
y0 = 1;

y = impEuler(a, t_h, y0, tsc);
current_best = y;

y_real = exp(a*t);
error(1) = norm(y_real - y);

%figure
%plot(t, y_real, t,y);
%legend('real sol', 'backward euler');

%%%%%%%%%%%% fine approximation in parallel %%%%%%%%%%%

%%%%%%%% PUT LOOP HERE %%%%%%
for iter=1:max_iter
    if(norm(y-y_real)<eps)
        fprintf('Converged after %d iterations', iter);
        break;
    end

% want to solve more exactly for each time interval
% first implement sequential with old values

% no of fine steps in each interval
no_fine_st = 10;
t_h_fine = T/(tsc*no_fine_st);

% think about overlap, last step should be first of new one
% so fine solution has tsc*(no_fine_st-1)+1 for last one entries
% but we actually want to keep track of all results at overlap

y_fine = impEuler(a, t_h_fine, (y(1:end-1))', no_fine_st);
% create corresponding x value matrix
x = zeros(tsc,no_fine_st+1);
for i=1:tsc
    x(i,:) = linspace((i-1)*t_h, i*t_h, no_fine_st+1);
end


%figure
%plot(t,y, x(1,:), y_fine(1,:), x(2,:), y_fine(2,:), x(3,:), y_fine(3,:), x(4,:), y_fine(4,:), x(5,:), y_fine(5,:));

% now we can see that we get jumps, compute jump
% need last element from one row and first element from second row
S = zeros(tsc+1,1);

for i=2:(tsc+1)
    S(i,1) = y_fine(i-1,end) - y(i);
end

% now also apply implicit Euler on jumps 

% change for more iterations
delta = zeros(1, tsc+1);

for i=2:tsc+1
    temp = impEuler(a, t_h, S(i-1,1)+delta(1,i-1), 1);
    delta(1,i) = temp(2);
end

% update coarse grid
y_old = y;
y(2:end) = y_fine(1:end,end)' + delta(1,2:end);

error(iter+1) = norm(y_real - y);

if(norm(y_real - current_best) > norm(y_real - y))
    current_best = y;
end

end

disp(error);

%%% fine grid Euler %%%
y_fine_approx = impEuler(a, t_h_fine, y(1), no_fine_st*tsc);

figure
plot(0:length(error)-1, error);
legend('error');
xlabel('iterations');

y_fine_vec = transpose(y_fine(:,1:end-1));
y_fine_vec = [y_fine_vec(:); y(end)];

t_fine = linspace(0,T,length(y_fine_vec));
y_real_fine = exp(a*t_fine);

figure
plot(t_fine, y_fine_approx, t, y, t_fine, y_fine_vec);
legend('implicit Euler fine', 'coarse Euler', 'parareal');

%figure
%plot(t, y, t, y_real, t_fine, y_fine_vec, t_fine, y_real_fine);
%legend('coarse approx', 'coarse real', 'fine approx', 'fine real');

%%%%%%%%%%%
