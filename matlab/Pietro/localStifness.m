function [A,L,D] = localStifness(p)
A=zeros(3);
L(:,1)= p(2,:)-p(1,:);
L(:,2)= p(3,:)-p(1,:);
Grad = [[-1,-1]',[1,0]',[0,1]'];

%G=[1.1278 0.8722;0.8722 1.1278]; % Ass 13

D=det(L);
area=0.5*abs(D);

for i = 1:3
    for j = 1:3
         A(i,j) = area*((L'\Grad(:,j)))'* (L'\Grad(:,i));
    end
end
end

