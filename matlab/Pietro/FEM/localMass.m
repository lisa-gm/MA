function [M] = localMass(p)
L(:,1)= p(2,:)-p(1,:);
L(:,2)= p(3,:)-p(1,:);
area=abs(det(L));
M=[2 1 1;1 2 1;1 1 2]*(1/24);
M=area*M;
end
