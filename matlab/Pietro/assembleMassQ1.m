%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% 2D MASS MATRIX PIETRO %%%%%%%%%%%%%%%%%%

function matrice=assembleMassQ1(M,N)

offset=M+1;

matrice=zeros((M+1)*(N+1));

%A=assembleReferenceMass(1);

A =  [ 1/9            1/18           1/18           1/36    
       1/18           1/9            1/36           1/18    
       1/18           1/36           1/9            1/18    
       1/36           1/18           1/18           1/9];

for J=1:N
    for I=1:M
        local2Global=[ (J-1)*offset+I (J-1)*offset+I+1 (J)*offset+I (J)*offset+I+1];
        
        for i=1:4
            for j=1:4
                matrice(local2Global(i),local2Global(j))=matrice(local2Global(i),local2Global(j))+A(i,j);
            end
        end
        
        
    end
end

end

