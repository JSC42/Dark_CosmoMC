function [r] = Integrate_E(z_dep_id, z_inj_id,m,dep_channel)
global T_Array
global Energy_Axis
global zp_Axis

% low-E axis
E1 = 100;
NE = 100;
E2 = Energy_Axis(1)-1;
le1 = log10(E1);
le2 = log10(E2);
dle = (le2 - le1)/(NE-1);
le = le1:dle:le2;
E_Low = 10.^le;
E_Low=E_Low';

E_Full = [E_Low;Energy_Axis];

L = dLdE(zp_Axis(z_dep_id)-1, m, E_Full);
T(:) = T_Array(1, dep_channel,z_dep_id, :, z_inj_id);
T=T';
% Apparantly low-E side can be extrapolated to 100 eV
T_Low=T(1)*ones(NE,1);
T_Full=[T_Low;T];
fun = T_Full.*L;

r = trapz(E_Full,fun);

end