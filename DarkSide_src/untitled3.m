cd /Users/cangtao/IHEPBox/Projects/GitHub/Dark_CosmoMC/DarkSide_src
clear

len = 2
LineWidth=2
PlotSize=18
id = 1000

clc
F1 = '/Users/cangtao/IHEPBox/Projects/GitHub/Dark_CosmoMC/HyRec/tmp_1.dat';
F2 = '/Users/cangtao/IHEPBox/Projects/GitHub/Dark_CosmoMC/HyRec/tmp_2.dat';
F1=load(F1);
F2=load(F2);
z1=F1(:,1);
x1=F1(:,2);
t1=F1(:,3);

z2=F2(:,1);
x2=F2(:,2);
t2=F2(:,3);

clf
subplot(1,2,1)
loglog(z1,x1,'k','LineWidth',LineWidth);hold on
loglog(z2,x2,'+k','LineWidth',LineWidth);hold on
xlabel('$z$','Interpreter','latex','FontSize',PlotSize,'FontName','Times');
ylabel('$x$','Interpreter','latex','FontSize',PlotSize,'FontName','Times');
set(gca,'FontSize',PlotSize,'Fontname','Times');

subplot(1,2,2)
loglog(z1,t1,'k','LineWidth',LineWidth);hold on
loglog(z2,t2,'+k','LineWidth',LineWidth);hold on
xlabel('$z$','Interpreter','latex','FontSize',PlotSize,'FontName','Times');
ylabel('$T$','Interpreter','latex','FontSize',PlotSize,'FontName','Times');
set(gca,'FontSize',PlotSize,'Fontname','Times');


G=6.67E-11;
Ms = 1.98847E30;
Ob = 0.048974681618697;
h=0.6766;
rhocr = 1.879E-26 * h^2;
km=1E3;
4*pi*G^2*Ms^2*Ob*rhocr/(1.44E6 * km^3)
