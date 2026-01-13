%% Function: plot 2D lines (in Tecplot ASCII) 
function TecplotLines(Data,time,Filename)
    [np,nvar] = size(Data);
    fid = fopen(Filename,'w');
    write_header = 'VARIABLES = time';
    for i = 1:nvar
       write_header = [write_header,',var',num2str(i)]; 
    end
    fprintf(fid,'%s',write_header);
    if length(time) == 1
        Data = [(0:time:(np-1)*time)',Data];
    elseif length(time) == np
        Data = [time', Data];
    else
        error('size of time does not matches sampling');
    end
    for i = 1:np
        fprintf(fid, '\n');
        for j = 1:nvar+1
            fprintf(fid, '%.4e ',Data(i,j));
        end
    end
    
    fclose(fid);
    fprintf('Successfully wrote file: %s\n', Filename);
    end

    