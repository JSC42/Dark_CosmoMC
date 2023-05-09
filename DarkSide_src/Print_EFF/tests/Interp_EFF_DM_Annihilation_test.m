cd /Users/cangtao/IHEPBox/Projects/GitHub/HyRec/DarkSide_src
clear

N = 100
reload = 1

clc
z_file = 'EFF/Redshift_1.txt';
m_file = 'EFF/M_ANN_1.txt';
EFF_file = 'EFF/Higgs_LyA_HMG.txt';
lm_axis=log10(load(m_file)/1E9);
lz_axis=log10(load(z_file));
EFF=load(EFF_file);
lm1 = min(lm_axis);
lm2 = max(lm_axis);
lz1 = min(lz_axis);
lz2 = max(lz_axis);

if reload
    for id = 1:N
        lm = lm1+rand*(lm2-lm1);
        lz = lz1+rand*(lz2-lz1);
        
        ! cp ~/FUNCTIONS/interp2_.m ./
        r1 = interp2_(lz_axis,lm_axis,EFF,lz,lm,1);
        ! rm interp2_.m
        M_GeV = 10^lm;
        z = 10^lz;
        
        OutFile = 'tmp.dat';
        delete(OutFile);
        FileID=fopen(OutFile,'a');
        fprintf(FileID,'%f\n',M_GeV);
        fprintf(FileID,'%f',z);
        fclose(FileID);
        ! ./a.out<tmp.dat>tmp_out.dat
        r2=load('tmp_out.dat');
        
        dif_1 = abs(1-r1/r2);
        dif_2 = abs(1-r2/r1);
        
        dif(id) = max(dif_1, dif_2);
        if dif(id) - dif(id) ~= 0
            dif(id) = 0;
        end
        
    end
    save tmp dif
else
    load tmp
end

max_dif = max(dif)
min_dif = min(dif)
avg_dif = sum(dif)/N
