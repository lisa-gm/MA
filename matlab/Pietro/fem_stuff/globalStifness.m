function [A,M]=globalStifness(file_name)
mesh=read_mesh(file_name);
Ne=mesh.number_of_elements; Nv=mesh.number_of_vertices;
vertices=mesh.vertices; elem=mesh.elements;
A=zeros(Nv);
M=zeros(Nv);

% elem gives the 3 vertices defining each triangle
for n=1:Ne
    % returns in both cases 3x3 matrix
    % compute elementwise, not basis function wise!
    % therefore run through all elements! 
    % and add it in the right spot
    Al=localStifness(vertices(:,elem(:,n))');
    Ml=localMass(vertices(:,elem(:,n))');
    
    for i=1:3
        for j=1:3
            % elem(i,n) gives a node
            A(elem(i,n),elem(j,n))=A(elem(i,n),elem(j,n))+Al(i,j);
            M(elem(i,n),elem(j,n))=M(elem(i,n),elem(j,n))+Ml(i,j);
        end
    end
end
end

% recipe: 
% go through all elements i.e. triangles
% for each triangle compute how the different vertices relate 
% to each other, using reference triangle
% hence only compute on that local triangle what int ( phi_i * phi_j) is
% if we look at triangle defined by vertices 1,5,7 
% we compute how they relate to each other 
% and then write into A(1,1) the integral, into A(1,5) etc. 
% but then vertex 5 might also be in an element of vertices 
% 5,9,17, then when we go through elements later, get there
% and then add additional part of integral
% so we dont ever assemble basis function in one but consider 
% everything elementwise!



