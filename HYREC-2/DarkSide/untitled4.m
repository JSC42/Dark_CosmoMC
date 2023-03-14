cd /Users/cangtao/IHEPBox/Projects/GitHub/HyRec/DarkSide_src

xfile = 'EFF/Energy_1.txt';
yfile = 'EFF/Redshift_1.txt';
zfile = 'EFF/Elec_HIon_HMG.txt';
x = load(xfile);
y = load(yfile);
z = load(zfile);

clf
h=pcolor(x,y,z);hold on
set(h, 'EdgeColor', 'none');
colormap('jet')
set(gca,'XScale','log') 
set(gca,'YScale','log') 
size(x)
size(y)
size(z)
