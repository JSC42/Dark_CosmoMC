function Printer(Input, Tab_Name)


OutFile = 'Tables.h';

Data=load(Input);
Size=size(Data);
if (Size(1) == 1) ||(Size(2) == 1)
    DIM = 1;
    % Allow for some redundency
    if ( (Data(1) < 0) || (Data(end)<0) )
        fprintf('Redundency not applipicable\n');
    end
    Data(1)=Data(1)*0.99;
    Data(end)=Data(end)*1.01;
    
else
    DIM = 2;
end

if DIM == 1
    nx=length(Data);
    FileID=fopen(OutFile,'a');
    
    fprintf(FileID,'double ');
    fprintf(FileID,Tab_Name);
    fprintf(FileID,'[');
    fprintf(FileID,'%i',nx);
    fprintf(FileID,'] = {');
    for xid=1:nx-1
        fprintf(FileID,'%E, ',Data(xid));
    end
    fprintf(FileID,'%E};\n',Data(nx));
    fprintf(FileID,'\n');
    fclose(FileID);
    
elseif DIM == 2
    ny = Size(1);
    nx = Size(2);
    FileID=fopen(OutFile,'a');
    % loger section
    fprintf(FileID,'// double ');
    fprintf(FileID,Tab_Name);
    fprintf(FileID,'[');
    fprintf(FileID,'%i]',ny);
    fprintf(FileID,'[');
    fprintf(FileID,'%i]\n',nx);
    
    fprintf(FileID,'double ');
    fprintf(FileID,Tab_Name);
    fprintf(FileID,'[');
    fprintf(FileID,'%i] = {\n',ny*nx);
    
    for yid=1:ny
        for xid=1:nx-1
            fprintf(FileID,'%E, ',Data(yid,xid));
        end
        if yid==ny
            fprintf(FileID,'%E};\n\n',Data(yid,nx));
        else
            fprintf(FileID,'%E,\n',Data(yid,nx));
        end
        
    end
    fprintf(FileID,'\n');
    
    fclose(FileID);
end
end