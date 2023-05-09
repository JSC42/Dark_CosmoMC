cd /Users/cangtao/cloud/GitHub/Dark_CosmoMC/camb
% An example for pretty 1D plot
clear

LineWidth=2
PlotSize=21

clc

F1 = 'LCDM_lensedCls.dat';
F2 = 'test_lensedCls.dat';
F3 = 'test_PBH_lensedCls.dat';

F1=importdata(F1)
F2=importdata(F2);
F3=importdata(F3);
F1=F1.data;
F2=F2.data;
F3=F3.data;
l1 = F1(:,1);
e1 = F1(:,3);

l2 = F2(:,1);
e2 = F2(:,3);

l3 = F3(:,1);
e3 = F3(:,3);

clf
loglog(l1,e1,'k','LineWidth',1.5*LineWidth);hold on
%loglog(l2,e2,'r','LineWidth',LineWidth);hold on
loglog(l3,e3,'b','LineWidth',LineWidth);hold on
xlabel('$\ell$','Interpreter','latex','FontSize',PlotSize,'FontName','Times');
%ylabel('$D^{\rm{EE}}_{\ell}$','Interpreter','latex','FontSize',PlotSize,'FontName','Times');
set(gca,'FontSize',PlotSize,'Fontname','Times');
LgD=legend('$\Lambda$CDM',...
    'Energy Injection');
set(LgD,'Interpreter','latex','Location','southeast','FontSize',PlotSize)
axis([2 2000 1e-3 1e2])
