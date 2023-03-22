cd /Users/cangtao/IHEPBox/Projects/GitHub/Dark_CosmoMC/HyRec
clear

LineWidth=2
PlotSize=20

clc

d0=load('tmp_1.dat');
zp0=d0(:,1) + 1;
x0=d0(:,2);
t0=d0(:,3);

d1=load('tmp_2.dat');
zp1=d1(:,1) + 1;
x1=d1(:,2);
t1=d1(:,3);

clf
subplot(1,2,1)
loglog(zp0,x0,'k','LineWidth',LineWidth);hold on
loglog(zp1,x1,'--r','LineWidth',LineWidth);hold on
xlabel('$z$','Interpreter','latex','FontSize',PlotSize,'FontName','Times');
ylabel('$x_e$','Interpreter','latex','FontSize',PlotSize,'FontName','Times');
set(gca,'FontSize',PlotSize,'Fontname','Times');
LgD=legend('LCDM',...
    'PBH');
set(LgD,'Interpreter','latex','Location','best','FontSize',PlotSize)

subplot(1,2,2)
loglog(zp0,t0,'k','LineWidth',LineWidth);hold on
loglog(zp1,t1,'--r','LineWidth',LineWidth);hold on
xlabel('$z$','Interpreter','latex','FontSize',PlotSize,'FontName','Times');
ylabel('$T_{\rm{m}}$','Interpreter','latex','FontSize',PlotSize,'FontName','Times');
set(gca,'FontSize',PlotSize,'Fontname','Times');
