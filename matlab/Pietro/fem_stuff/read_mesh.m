function mesh = read_mesh( file_name )


    file_id = fopen(file_name);
    data_line = fgetl(file_id);
    
    data_mesh = sscanf(data_line,'%i');
    
    mesh.number_of_vertices = data_mesh(1);
    mesh.number_of_elements = data_mesh(2);
    mesh.number_of_boundary_sides = data_mesh(3);
    
    
    
    mesh.vertices=[];
    
    for i=1:mesh.number_of_vertices
        data_line = fgetl(file_id);
        variable = sscanf(data_line,'%f');
        
        mesh.vertices=[mesh.vertices variable(1:2)];
        
        mesh.vertices_flag(i) = int16(variable(3));
        
    end

    mesh.elements=[];
    
    for i=1:mesh.number_of_elements
        data_line = fgetl(file_id);
        variable = sscanf(data_line,'%i');
        
        mesh.elements=[mesh.elements variable(1:3)];
        
        mesh.elements_flag(i) = variable(4);
        
    end

    mesh.boundary_sides=[];
    
    for i=1 : mesh.number_of_boundary_sides
        
        data_line = fgetl(file_id);
        variable = sscanf(data_line,'%i');
        
        mesh.boundary_sides=[mesh.boundary_sides variable(1:2)];
        mesh.boundary_sides_flag(i) = variable(3);

        
    end

    fclose(file_id);
    