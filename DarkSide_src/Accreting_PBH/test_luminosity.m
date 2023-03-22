cd /Users/cangtao/IHEPBox/Projects/GitHub/Dark_CosmoMC/DarkSide_src/Accreting_PBH/
clear

reload = 1
N=100

clc
Initialise
lm = log10(Mass_Axis);
lz = log10(zp_Axis-1);
lz1 = min(lz);
lz2 = max(lz);
lm1 = min(lm);
lm2 = max(lm);
len_m = lm2-lm1;
len_z = lz2-lz1;

if reload
    tic
    for id = 1:N
        lm_ = lm1+rand*len_m;
        
        %lz_ = lz1+rand*len_z;
        lz_ = lz(10);
        
        m=10^lm_;
        z=10^lz_;
        r1 = Luminosity(z,m);
        
        cd ../../HyRec/
        File = 'tmp_in.dat';
        delete(File)
        FileID=fopen(File,'a');
        fprintf(FileID,'%E\n',m);
        fprintf(FileID,'%E',z);
        fclose(FileID);
        ! ./a.out<tmp_in.dat>tmp_out.dat
        r2 = load('tmp_out.dat');
        cd ../DarkSide_src/Accreting_PBH/
        dif(id) = (r1-r2)/r1;
    end
    toc
    save tmp dif
else
    load tmp
end

abs_dif = abs(dif);
max(abs_dif)
avg_dif = sum(dif)/N
