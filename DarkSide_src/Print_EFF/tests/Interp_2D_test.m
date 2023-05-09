cd /Users/cangtao/IHEPBox/Projects/GitHub/HyRec/DarkSide_src

reload = 1
N = 300

clc
z_file = 'EFF/Redshift_1.txt';
Ek_file = 'EFF/Energy_1.txt';
EFF_file = 'EFF/Phot_LyA_HMG.txt';

lek_axis=log10(load(Ek_file));
lzp_axis=log10(load(z_file));
EFF=load(EFF_file);
le1 = min(lek_axis);
le2 = max(lek_axis);
lz1 = min(lzp_axis);
lz2 = max(lzp_axis);

if reload
    for id = 1:N
        le = le1+rand*(le2-le1);
        lzp = lz1+rand*(lz2-lz1);
        
        ! cp ~/FUNCTIONS/interp2_.m ./
        r1 = interp2_(lzp_axis,lek_axis,EFF,lzp,le,1);
        ! rm interp2_.m
        Ek = 10^le;
        zp = 10^lzp;
        
        OutFile = 'tmp.dat';
        delete(OutFile);
        FileID=fopen(OutFile,'a');
        fprintf(FileID,'%f\n',Ek);
        fprintf(FileID,'%f',zp);
        fclose(FileID);
        ! ./a.out<tmp.dat>tmp_out.dat
        r2=load('tmp_out.dat');
        
        dif_1 = abs(1-r1/r2);
        dif_2 = abs(1-r2/r1);
        
        dif(id) = max(dif_1, dif_2);
    end
    save tmp dif
else
    load tmp
end

max_dif = max(dif)
min_dif = min(dif)
avg_dif = sum(dif)/N
