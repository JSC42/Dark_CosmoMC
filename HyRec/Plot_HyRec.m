cd /Users/cangtao/cloud/GitHub/Dark_CosmoMC/HyRec
clear

LineWidth=2
PlotSize=20

clc

d0=load('tmp_out.dat');
zp0=d0(:,1) + 1;
x0=d0(:,2);
t0=d0(:,3);

clf
subplot(1,2,1)
loglog(zp0,x0,'k','LineWidth',LineWidth);hold on
xlabel('$z$','Interpreter','latex','FontSize',PlotSize,'FontName','Times');
ylabel('$x_e$','Interpreter','latex','FontSize',PlotSize,'FontName','Times');
set(gca,'FontSize',PlotSize,'Fontname','Times');
%axis([0 2000 1e-4 1.4])

subplot(1,2,2)
loglog(zp0,t0,'k','LineWidth',LineWidth);hold on
xlabel('$z$','Interpreter','latex','FontSize',PlotSize,'FontName','Times');
ylabel('$T_{\rm{k}}$','Interpreter','latex','FontSize',PlotSize,'FontName','Times');
set(gca,'FontSize',PlotSize,'Fontname','Times');
%axis([0 2000 1e1 1e4])
