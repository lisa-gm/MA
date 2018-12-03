% change implementation so that it can run in parallel
% then y0 is a vector, each entry corresponds to a separate problem

function y = impEuler(A, h, y0, steps)

% y becomes matrix
y = zeros(size(y0,1),steps+1);
y(:,1) = y0;

I = eye(size(y0,1));

for i=2:(steps+1)
    % we want: y(i) = y(i-1) + h * f(i,y(i))
    % assume f = A*y, then y(i) = y(i-1) + h**y(i)
    % solve for y(i):
    y(:,i) = A * y(:,i-1);
end

end