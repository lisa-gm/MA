function u = heat_1D_anal( u0,x,t,N)

% u0 : initial (constant) value
% n : bigger n better solution 

u = 0;

for n = 1:2:2*N   
    u = u + (4*u0/pi) * (sin(n*pi*x)*(exp(-n^2*pi^2*t))/n);    
end

end

