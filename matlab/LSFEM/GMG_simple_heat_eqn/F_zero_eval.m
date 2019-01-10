% **************************************************** %
% ************ DEFINING A ZERO FORCING TERM ********** %
% **************************************************** %

% returns zeros of the correct size, for heat eqn
function [f,df,ddf] = F_zero_eval(u)
   f = 0*ones(size(u)); 
   df = 0*ones(size(u));  
   ddf = 0*ones(size(u));
end