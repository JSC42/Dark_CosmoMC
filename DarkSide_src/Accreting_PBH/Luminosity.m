function [r] = Luminosity(z,m)
% PBH luminosity in eV/s
% Speed: 4/10^5 s per call

global Luminosity_Table
global Mass_Axis
global zp_Axis
nz = length(zp_Axis);
nm = length(Mass_Axis);
zp = 1+z;
r=Interp_2D(zp_Axis,Mass_Axis,1+z,m,Luminosity_Table,2);

end