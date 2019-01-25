% **************************************************** %
% ************* DEFINING THE FORCING TERM ************ %
% **************************************************** %

function [f,df,ddf] = F_eval(u)
%    f = 0*ones(size(u)); %u.*u; 
%    df = 0*ones(size(u)); % 2.*u; %
%    ddf = 0*ones(size(u));
  alpha = 1; 
  th = 0.1;
  f = alpha*(-u.*(u-th).*(u-1));
  df = alpha*(-((u-th).*(u-1) + u.*(u-1) + u.*(u-th)));
  ddf = alpha*(-((u-th)+(u-1) + u+(u-1) + u+(u-th)));

end