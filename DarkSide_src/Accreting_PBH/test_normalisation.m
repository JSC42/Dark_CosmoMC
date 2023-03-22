cd /Users/cangtao/IHEPBox/Projects/GitHub/Dark_CosmoMC/DarkSide_src/Accreting_PBH/
clear

z = 102
m = 1E2
le1 = 1
le2 = 12.9
dle = 0.01

clc
le = le1:dle:le2;
E=10.^le;
L = dLdE(z,m,E);
r1 = trapz(E,L)
r2 = Luminosity(z,m)

Initialise
tic
Get_EFF(6,1E5,1)
toc
