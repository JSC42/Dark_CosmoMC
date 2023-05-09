function Print_EFF(dep_channel)

Path = '../Print_EFF/EFF/';
Suffix = '.txt';
Head = 'PBH_Accretion_Naked_';
if dep_channel==1
    Channel = 'HIon';
elseif dep_channel==3
    Channel = 'LyA';
elseif dep_channel==4
    Channel = 'Heat';
else
    error('Wrong channel setting.\n')
end

S = string({Path,Head,Channel,Suffix});
FileName = char(join(S,''))
delete(FileName)

load ./data/EFF.mat
Data(:,:) = EFF(dep_channel, :,:);
N=size(Data);
nz = N(1);
nm = N(2);

FileID=fopen(FileName,'w');
for zid = 1:nz
    for mid = 1:nm
        f = Data(zid,mid);
        fprintf(FileID,'%E	',f);
    end
    fprintf(FileID,'\n');
end
fclose(FileID);

end