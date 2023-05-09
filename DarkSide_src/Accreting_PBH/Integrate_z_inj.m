function [r]=Integrate_z_inj(z_dep_id,m,dep_channel)
global zp_Axis
nz = length(zp_Axis);

for z_inj_id = 1:nz
    LevelE(z_inj_id,1) = Integrate_E(z_dep_id, z_inj_id,m,dep_channel);    
end

z = zp_Axis-1;
H=Hubble(z);
fun = LevelE./(H.*(1+z));
r = trapz(z,fun);

end