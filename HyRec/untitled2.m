cd /Users/cangtao/cloud/GitHub/Dark_CosmoMC/HyRec
clear

LineWidth=4
PlotSize=30

clc

d0=load('tmp_1.dat');
zp0=d0(:,1) + 1;
x0=d0(:,2);
t0=d0(:,3);

d1=load('tmp_2.dat');
zp1=d1(:,1) + 1;
x1=d1(:,2);
t1=d1(:,3);

d2=load('tmp_3.dat');
zp2=d2(:,1) + 1;
x2=d2(:,2);
t2=d2(:,3);

clf
%subplot(2,1,1)
subplot(1,2,1)
loglog(zp0,x0,'k','LineWidth',LineWidth);hold on
%loglog(zp1,x1,'r','LineWidth',LineWidth);hold on
loglog(zp2,x2,'b','LineWidth',LineWidth);hold on
xlabel('$z$','Interpreter','latex','FontSize',PlotSize,'FontName','Times');
%ylabel('Ionisation Fraction','Interpreter','latex','FontSize',PlotSize,'FontName','Times');
set(gca,'FontSize',PlotSize,'Fontname','Times');
title('Ionisation Fraction','Interpreter','latex','FontSize',PlotSize)
LgD=legend('$\Lambda$CDM',...
    'Energy Injection');
set(LgD,'Interpreter','latex','Location','northwest','FontSize',PlotSize)
axis([12 2000 1e-4 1.4])

%subplot(2,1,2)
subplot(1,2,2)
loglog(zp0,t0,'k','LineWidth',LineWidth);hold on
%loglog(zp1,t1,'r','LineWidth',LineWidth);hold on
loglog(zp2,t2,'b','LineWidth',LineWidth);hold on
xlabel('$z$','Interpreter','latex','FontSize',PlotSize,'FontName','Times');
% ylabel('Gas Temperature (k)','Interpreter','latex','FontSize',PlotSize,'FontName','Times');
title('Gas Temperature (k)','Interpreter','latex','FontSize',PlotSize)
set(gca,'FontSize',PlotSize,'Fontname','Times');
axis([12 2000 1e1 1e4])

