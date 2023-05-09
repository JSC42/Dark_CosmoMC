cd /Users/cangtao/IHEPBox/Projects/GitHub/HyRec/DarkSide_src
clear

reload = 0

clc
EFF_file = 'EFF/Phot_HIon_HMG.txt';
f=load(EFF_file);
Size=size(f);
nz=Size(1);
nm=Size(2);

if reload
    for zid = 1:nz
        for mid = 1:nm
            OutFile = 'tmp.dat';
            delete(OutFile);
            FileID=fopen(OutFile,'a');
            fprintf(FileID,'%i\n',mid-1);
            fprintf(FileID,'%i',zid-1);
            fclose(FileID);
            ! ./a.out<tmp.dat>tmp_out.dat
            r1 = f(zid,mid);
            r2 = load('tmp_out.dat');
            dif_1 = abs(1-r1/r2);
            dif_2 = abs(1-r2/r1);
            
            dif(zid,mid) = max(dif_1, dif_2);
            
        end
        status = zid/nz
    end
    save tmp dif
else
    load tmp
end

max_dif=max(max(dif))
min_dif=min(min(dif))

avg_dif = sum(sum(dif))/(nz*nm)
