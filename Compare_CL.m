cd /Users/cangtao/IHEPBox/Projects/GitHub/Dark_CosmoMC/

id = 4
LineWidth=2
PlotSize=18

clc
f1='camb/tmp_0.dat';
f2='camb/tmp_1.dat';

d1=importdata(f1);
d2=importdata(f2);
d1=d1.data();
d2=d2.data();

l1=d1(:,1);
l2=d2(:,1);
x1=d1(:,id);
x2=d2(:,id);

clf
loglog(l1,x1,'k','LineWidth',LineWidth);hold on
loglog(l2,x2,'r','LineWidth',LineWidth);hold on
xlabel('$x$','Interpreter','latex','FontSize',PlotSize,'FontName','Times');
ylabel('$y$','Interpreter','latex','FontSize',PlotSize,'FontName','Times');
set(gca,'FontSize',PlotSize,'Fontname','Times');

