function [A,M]=globalStifness(file_name)
mesh=read_mesh(file_name);
Ne=mesh.number_of_elements; Nv=mesh.number_of_vertices;
vertices=mesh.vertices; elem=mesh.elements;
elem(:,1)
A=zeros(Nv);
M=zeros(Nv);

for n=1:Ne
    % returns in both cases 3x3 matrix
    Al=localStifness(vertices(:,elem(:,n))');
    Ml=localMass(vertices(:,elem(:,n))');
    
    for i=1:3
        for j=1:3
            % write into global matrix
            % where do we write what? 
            % 
            A(elem(i,n),elem(j,n))=A(elem(i,n),elem(j,n))+Al(i,j);
            M(elem(i,n),elem(j,n))=M(elem(i,n),elem(j,n))+Ml(i,j);
        end
    end
end
end



