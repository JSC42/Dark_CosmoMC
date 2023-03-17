cd /Users/cangtao/IHEPBox/Projects/GitHub/HyRec/DarkSide_src

N = 100

clc
z_file = 'EFF/Redshift_1.txt';
m_file = 'EFF/Energy_1.txt';
EFF_file = 'EFF/Phot_HIon_HMG.txt';
m_axis=log10(load(m_file)/1E9);
z_axis=log10(load(z_file));
EFF=load(EFF_file);
m1 = min(m_axis);
m2 = max(m_axis);
z1 = min(z_axis);
z2 = max(z_axis);

lm = m1+rand*(m2-m1);
lz = z1+rand*(z2-z1);

for kk=1:N
Seed1 = rand;
Seed2 = rand;
mid = 1+round(length(m_axis)*Seed1);
zid = 1+round(length(z_axis)*Seed2);
mid=min(length(m_axis),mid);
zid=min(length(z_axis),zid);
mid=max(1,mid);
zid=max(1,zid);
lm = m_axis(mid);
lz = z_axis(zid);
r0 = EFF(zid,mid);
! cp ~/FUNCTIONS/interp2_.m ./
r1 = interp2_(z_axis,m_axis,EFF,lz,lm,1);
! rm interp2_.m
dif = abs(1-r0/r1);
if dif > 1E-4
    ErrorVar
end
end
