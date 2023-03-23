cd /Users/cangtao/IHEPBox/Projects/GitHub/Dark_CosmoMC/DarkSide_src/Accreting_PBH/
clear

reload = 0
nm = 400

clc
Initialise
nz = length(zp_Axis);
lm1 = log10(Mass_Axis(1));
lm2 = log10(Mass_Axis(end));
dlm = (lm2 - lm1)/(nm-1);
lm = lm1:dlm:lm2;
M_Axis = 10.^lm;

% ---- Calculate EFF ----
if reload
    tic
    for mid = 1:nm
        for zid = 1:nz
            m = M_Axis(mid);
            for cid = 1:5
                % zid = nz tends to give interpolation error, just use nearest nz-1
                if zid == nz
                    r = Get_EFF(nz-1, m, cid);
                else
                    r = Get_EFF(zid, m, cid);
                end
                EFF(cid,zid,mid) = r;
            end
        end
        status = mid/nm
    end
    toc
    save ./data/EFF.mat EFF M_Axis zp_Axis
else
    load ./data/EFF.mat
end

% ---- Print Mass Axis ----
Mass_Axis_File = '../Print_EFF/EFF/PBH_Accretion_Mass_Axis.txt';
delete(Mass_Axis_File)
FileID=fopen(Mass_Axis_File,'a');
for mid = 1:nm
    fprintf(FileID,'%E\n',M_Axis(mid));
end
fclose(FileID);

% ---- Print EFF ----
Print_EFF(1);
Print_EFF(3);
Print_EFF(4);
