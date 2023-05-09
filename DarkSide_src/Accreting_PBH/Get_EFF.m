function [r]=Get_EFF(z_dep_id, m, dep_channel)
global zp_Axis
x = log(zp_Axis);
x1 = x(1:end-1);
x2 = x(2:end);
dx = sum(abs(x2-x1))/length(x2);

z = zp_Axis(z_dep_id);
H = Hubble(z);
Lb = Luminosity(z,m);
Level_z_inj=Integrate_z_inj(z_dep_id,m,dep_channel);

r = H*Level_z_inj/Lb/dx;

end