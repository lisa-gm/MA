function [A,M]=globalStifness(file_name)
mesh=read_mesh(file_name);
Ne=mesh.number_of_elements; Nv=mesh.number_of_vertices;
vertices=mesh.vertices; elem=mesh.elements;
A=zeros(Nv);
M=zeros(Nv);

for n=1:Ne
    Al=localStifness(vertices(:,elem(:,n))');
    Ml=localMass(vertices(:,elem(:,n))');
    
    for i=1:3
        for j=1:3
            A(elem(i,n),elem(j,n))=A(elem(i,n),elem(j,n))+Al(i,j);
            M(elem(i,n),elem(j,n))=M(elem(i,n),elem(j,n))+Ml(i,j);
        end
    end
end
end



