global Luminosity_Table
global Mass_Axis
global zp_Axis
global Energy_Axis
global T_Array
Luminosity_Table = load('data/Luminosity_Table.txt');
Mass_Axis=load('data/Mbh_Axis.txt');
zp_Axis=load('../Print_EFF/EFF/Redshift_1.txt');
Energy_Axis=load('../Print_EFF/EFF/Energy_1.txt');
T_Array = h5read('data/Tracy.h5','/T_Array');

% f='data/Tracy.h5'
% h5disp(f)
% 
% x=log10(zp_Axis);
% x1=x(1:end-1);
% x2=x(2:end);
% dx = x2-x1;
