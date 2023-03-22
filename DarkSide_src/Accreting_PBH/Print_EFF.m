function Print_EFF(Particle_Channel, dep_channel, process)

Path = '/Users/cangtao/IHEPBox/Projects/GitHub/Dark_CosmoMC/DarkSide_src/Print_EFF/EFF/';
suffix = '.txt';

Particle_ID = Match_Channel(Particle_Channel);

if process==1
    Process_Name = 'DEC';
elseif process==2
    Process_Name = 'HMG';
end

if dep_channel==1
    Dep_Channel_Name = 'HIon_';
elseif dep_channel==3
    Dep_Channel_Name = 'LyA_';
elseif dep_channel==4
    Dep_Channel_Name = 'Heat_';
else
    error('wrong dep_channel.\n');
end

if Particle_Channel == 3
    Particle_Name = 'Higgs_';
elseif Particle_Channel == 4
    Particle_Name = 'Muon_';
elseif Particle_Channel == 5
    Particle_Name = 'Tau_';
elseif Particle_Channel == 6
    Particle_Name = 'Q_';
elseif Particle_Channel == 7
    Particle_Name = 'Charm_';
elseif Particle_Channel == 8
    Particle_Name = 'Bottom_';
elseif Particle_Channel == 9
    Particle_Name = 'Top_';
elseif Particle_Channel == 10
    Particle_Name = 'W_';
elseif Particle_Channel == 11
    Particle_Name = 'Z_';
elseif Particle_Channel == 12
    Particle_Name = 'Gluon_';
else
    error('Wrong Particle Channel.\n');
end

S = string({Path,Particle_Name,Dep_Channel_Name,Process_Name,suffix});
FileName = char(join(S,''));

% clean up
delete(FileName)
load ./Data/RESULT.mat
% this is the data to print
fc(:,:) = EFF(Particle_ID,dep_channel,:,process,:);

N = size(fc);
nz = N(2);
nm = N(1);

FileID=fopen(FileName,'a');
for zid = 1:nz
    for mid = 1:nm
        f = fc(mid,zid);
        fprintf(FileID,'%E	',f);
    end
    fprintf(FileID,'\n');
end
fclose(FileID);
end