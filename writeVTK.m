function writeVTK(filename, vertices, faces, varargin)
    % Open file for writing
    fileID = fopen(filename, 'w');
    
    % Write VTK header
    fprintf(fileID, '# vtk DataFile Version 3.0\n');
    fprintf(fileID, 'PLY to VTK conversion\n');
    fprintf(fileID, 'ASCII\n');
    fprintf(fileID, 'DATASET POLYDATA\n');
    
    % Write vertices
    fprintf(fileID, 'POINTS %d float\n', size(vertices, 1));
    fprintf(fileID, '%f %f %f\n', vertices');
    
    % Write faces
    fprintf(fileID, 'POLYGONS %d %d\n', size(faces, 1), size(faces, 1) * 4);
    fprintf(fileID, '%d %d %d %d\n', [3 * ones(size(faces, 1), 1), faces - 1]'); % VTK uses 0-based indexing
    
    % Check if scalar or vector fields are provided
    if ~isempty(varargin)
        fprintf(fileID, 'POINT_DATA %d\n', size(vertices, 1));
    end
    
    % Process scalar and vector fields
    i = 1;
    while i <= length(varargin)
        % Read the name and data
        fieldName = varargin{i};
        fieldData = varargin{i+1};
        
        % Check if the data is a scalar or vector field
        if size(fieldData, 2) == 1
            % Scalar field
            fprintf(fileID, 'SCALARS %s float 1\n', fieldName);
            fprintf(fileID, 'LOOKUP_TABLE default\n');
            fprintf(fileID, '%f\n', fieldData);
        elseif size(fieldData, 2) == 3
            % Vector field
            fprintf(fileID, 'VECTORS %s float\n', fieldName);
            fprintf(fileID, '%f %f %f\n', fieldData');
        else
            error('Field "%s" must be either a scalar (Nx1) or vector (Nx3).', fieldName);
        end
        
        % Move to the next pair of arguments
        i = i + 2;
    end
    
    % Close file
    fclose(fileID);
end
