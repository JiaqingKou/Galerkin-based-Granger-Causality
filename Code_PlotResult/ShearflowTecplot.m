%% Function: visualiztion of the shear flow model (in Tecplot ASCII format)
%   INPUT
%     Vec9：coefficients of 9 modes in vector, size(1,9)
%     Dt  : time step of snapshots
%     Filename: output file name, i.e. 'flow'
% 
%   OUTPUT
%     plot 3D fluid snapshots in Tecplot ASCII format
%
function ShearflowTecplot(Vec9, Dt, Filename)
    [snap, ~] = size(Vec9);
    gridx = 10;         % number of grids in x direction
    gridy = 20;         % number of grids in y direction
    gridz = 20;         % number of grids in z direction
    % Tecplot dimension mapping: I = ny, J = nx, K = nz
    I = gridy + 1;      % number of nodes in y direction -> Tecplot I
    J = gridx + 1;      % number of nodes in x direction -> Tecplot J
    K = gridz + 1;      % number of nodes in z direction -> Tecplot K
    % Tecplot file header
    file_header1 = 'TITLE = "Visualization of the 9-mode shear flow"';
    file_header2 = 'VARIABLES = "x","y","z","u","v","w","velocity"';
    file_header4 = ['I=', num2str(I), ', J=', num2str(J), ', K=', num2str(K), ', ZONETYPE=Ordered'];
   
    disp('# Tecplot File Generating...0%');
    for i = 1:snap
        if snap > 1
            % append frame number if multiple snapshots
            file_name = [Filename, '_', num2str(i), '.dat'];
        else
            % single snapshot
            file_name = [Filename, '.dat'];  
        end
        
        file_header3 = ['ZONE STRANDID=', num2str(i), ', SOLUTIONTIME=', num2str(i*Dt), ', DATAPACKING=POINT'];
        
        % get velocity field data
        [X, Y, Z, u, v, w] = ShearflowVelocity(Vec9(i, :), gridx, gridy, gridz);
        vel = sqrt(u.^2 + v.^2 + w.^2);     % velocity magnitude
        
        % reshape data to Tecplot required order: (I,J,K) -> innermost loop is I
        X = permute(X, [2, 1, 3]);  % Adjust dimension order: (x,y,z) -> (y,x,z)
        Y = permute(Y, [2, 1, 3]);  % y becomes first dimension (I)
        Z = permute(Z, [2, 1, 3]);  % x becomes second dimension (J)
        u = permute(u, [2, 1, 3]);
        v = permute(v, [2, 1, 3]);
        w = permute(w, [2, 1, 3]);
        vel = permute(vel, [2, 1, 3]);
        
        % flatten to 1D vector (in Tecplot order)
        X = X(:);
        Y = Y(:);
        Z = Z(:);
        u = u(:);
        v = v(:);
        w = w(:);
        vel = vel(:);

        % open file for writing
        fid = fopen(file_name, 'w');
        if fid == -1
            error(['Cannot open file: ', file_name]);
        end
        
        % write headers
        fprintf(fid, '%s\n%s\n%s\n%s\n', file_header1, file_header2, file_header3, file_header4);
        
        % write data points
        for j = 1:length(X)
            fprintf(fid, '%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\n', ...
                X(j), Y(j), Z(j), u(j), v(j), w(j), vel(j));
        end
        
        fclose(fid);
        
        % display progress
        if mod(i, 10) == 0
            fprintf('# Tecplot File Generating... %.1f%%\n', i/snap*100);
        end
    end
    disp('# Tecplot File Generated !');
end

