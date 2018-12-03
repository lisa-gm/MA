% **************************************************** %
% ************* DEFINING THE FORCING TERM ************ %
% **************************************************** %

function [f,df,ddf] = F_eval(u)
%   f = u.*u; %-3*ones(size(u));
%   df = 2.*u; %1*ones(size(u));
%   ddf = 2*ones(size(u));
  alpha = 1; 
  th = 0.1;
  f = alpha*(-u.*(u-th).*(u-1));
  df = alpha*(-((u-th).*(u-1) + u.*(u-1) + u.*(u-th)));
  ddf = alpha*(-((u-th)+(u-1) + u+(u-1) + u+(u-th)));

end