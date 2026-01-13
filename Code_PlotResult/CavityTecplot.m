%% Function: visualiztion of the cavity flow (in Tecplot ASCII)
%   Write structured grid data to Tecplot format DAT file
%   Supports 2D structured grid for contour plotting
%   
%   Input parameters:
%       FILENAME: Output filename (without extension)
%       DATACELL: Cell array containing variable names and data in pairs
%                 Format: {'VarName1', data1, 'VarName2', data2, ...}
%       NI: Number of points in I-direction (x-direction)
%       NJ: Number of points in J-direction (y-direction)
%   
%   Example usage:
%       dataCell = {'X', gridx, 'Y', gridy, 'Ux', ux, 'Uy', uy, 'Velocity', velocity};
%       CavityTecplot('output', dataCell, 129, 129);
%
function CavityTecplot(filename, dataCell, ni, nj)
    % Check input arguments
    if nargin < 4
        error('Filename, data cell array, NI, and NJ must be provided');
    end
    
    if ~iscell(dataCell)
        error('Second input must be a cell array');
    end
    
    if mod(length(dataCell), 2) ~= 0
        error('Data cell must contain variable names and data in pairs');
    end
    
    % Validate grid dimensions
    totalPoints = ni * nj;
    dataLength = length(dataCell{2});
    
    if dataLength ~= totalPoints
        error('Data length (%d) does not match grid dimensions %d x %d = %d', ...
              dataLength, ni, nj, totalPoints);
    end
    
    % Extract variable names and data
    numVars = length(dataCell) / 2;
    varNames = cell(1, numVars);
    data = cell(1, numVars);
    
    for i = 1:numVars
        nameIndex = (i-1)*2 + 1;
        dataIndex = nameIndex + 1;
        
        if ~ischar(dataCell{nameIndex})
            error('Variable name at position %d must be a string', nameIndex);
        end
        
        currentData = dataCell{dataIndex};
        if length(currentData) ~= totalPoints
            error('Data vector "%s" has incorrect length', dataCell{nameIndex});
        end
        
        varNames{i} = dataCell{nameIndex};
        data{i} = currentData(:); % Ensure column vector
    end
    
    % Create full filename
    full_filename = [filename, '.dat'];
    fid = fopen(full_filename, 'wt');
    
    if fid == -1
        error('Cannot create file: %s', full_filename);
    end
    
    try
        % Write Tecplot header for structured grid
        fprintf(fid, 'TITLE = "%s"\n', filename);
        fprintf(fid, 'VARIABLES = ');
        
        % Write variable names
        for i = 1:numVars
            if i == numVars
                fprintf(fid, '"%s"', varNames{i});
            else
                fprintf(fid, '"%s", ', varNames{i});
            end
        end
        fprintf(fid, '\n');
        
        % Write structured grid zone information
        fprintf(fid, 'ZONE T="Structured Grid", I=%d, J=%d, F=POINT\n', ni, nj);
        
        % Write data in structured order
        for j = 1:nj
            for i = 1:ni
                % Calculate index in the flattened array
                idx = (j-1)*ni + i;
                
                for k = 1:numVars
                    value = data{k}(idx);
                    if k == numVars
                        fprintf(fid, '%14.6e\n', value);
                    else
                        fprintf(fid, '%14.6e ', value);
                    end
                end
            end
        end
        
        fclose(fid);
        fprintf('Successfully wrote %dx%d structured grid with %d variables to: %s\n', ...
                ni, nj, numVars, full_filename);
        
    catch ME
        if fid ~= -1
            fclose(fid);
        end
        rethrow(ME);
    end
end