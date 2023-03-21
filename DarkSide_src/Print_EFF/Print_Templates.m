cd /Users/cangtao/IHEPBox/Projects/GitHub/Dark_CosmoMC/DarkSide_src/Print_EFF
clear

clc
F = 'EFF/LCDM_Xe_Tm_Template.dat';
F=load(F);
z=F(:,1);
x=F(:,2);
t=F(:,3);
z2=load('EFF/Redshift_2.txt')-1;

for zid = 1:length(z2)
    x2(zid) = interp1(z,x,z2(zid));
    t2(zid) = interp1(z,t,z2(zid));
end


F_xe = 'EFF/Xe_Template.txt';
F_T = 'EFF/Tm_Template.txt';
delete(F_xe)
delete(F_T)

Fx=fopen(F_xe,'a');
Ft=fopen(F_T,'a');
zp = 1+z2;
for zid=1:length(z2)
    fprintf(Fx,'%f\n',x2(zid));
    fprintf(Ft,'%f\n',t2(zid));
end

fclose(Fx);
fclose(Ft);
