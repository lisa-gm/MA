classdef solution_xv < solution
    % 
    % Solution of a system with position and velocity component
    %
    % x' = f1(x,v)
    % v' = f2(x,v)
    %
      
    methods(Abstract=true)
       
        %
        f1(obj)
        
        %
        f2(obj)
        
        
        
    end
end