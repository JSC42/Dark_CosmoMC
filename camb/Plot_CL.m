cd /Users/cangtao/cloud/GitHub/Dark_CosmoMC/camb
% An example for pretty 1D plot
clear

LineWidth=2
PlotSize=25

clc

F1 = 'LCDM_lensedCls.dat';
F2 = 'test_lensedCls.dat';

F1=importdata(F1)
F2=importdata(F2);
F1=F1.data;
F2=F2.data;
l1 = F1(:,1);
t1 = F1(:,2);
e1 = F1(:,3);
x1 = F1(:,5);

l2 = F2(:,1);
t2 = F2(:,2);
e2 = F2(:,3);
x2 = F2(:,5);

clf
subplot(1,3,1)
semilogx(l1,t1,'k','LineWidth',LineWidth);hold on
semilogx(l2,t2,'r','LineWidth',LineWidth);hold on
xlabel('$\ell$','Interpreter','latex','FontSize',PlotSize,'FontName','Times');
ylabel('$D^{\rm{TT}}_{\ell}$','Interpreter','latex','FontSize',PlotSize,'FontName','Times');
set(gca,'FontSize',PlotSize,'Fontname','Times');
LgD=legend('LCDM',...
    'DM');
set(LgD,'Interpreter','latex','Location','best','FontSize',PlotSize)

subplot(1,3,2)
semilogx(l1,x1,'k','LineWidth',LineWidth);hold on
semilogx(l2,x2,'r','LineWidth',LineWidth);hold on
xlabel('$\ell$','Interpreter','latex','FontSize',PlotSize,'FontName','Times');
ylabel('$D^{\rm{TE}}_{\ell}$','Interpreter','latex','FontSize',PlotSize,'FontName','Times');
set(gca,'FontSize',PlotSize,'Fontname','Times');

subplot(1,3,3)
semilogx(l1,e1,'k','LineWidth',LineWidth);hold on
semilogx(l2,e2,'r','LineWidth',LineWidth);hold on
xlabel('$\ell$','Interpreter','latex','FontSize',PlotSize,'FontName','Times');
ylabel('$D^{\rm{EE}}_{\ell}$','Interpreter','latex','FontSize',PlotSize,'FontName','Times');
set(gca,'FontSize',PlotSize,'Fontname','Times');
