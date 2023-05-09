cd ~/cloud/GitHub/Dark_CosmoMC/DarkSide_src/Accreting_PBH/
clear

LineWidth=2
PlotSize=20

cid = 4
mid = 201

clc
load data/EFF.mat
f1(:) = EFF(1,:,mid);
f3(:) = EFF(3,:,mid);
f4(:) = EFF(4,:,mid);
m=M_Axis(mid)

% -- Get plot
clf
loglog(zp_Axis,f1,'k','LineWidth',LineWidth);hold on
loglog(zp_Axis,f3,'r','LineWidth',LineWidth);hold on
loglog(zp_Axis,f4,'b','LineWidth',LineWidth);hold on
xlabel('$1+z$','Interpreter','latex','FontSize',PlotSize,'FontName','Times');
ylabel('$f_{\rm{c}}$','Interpreter','latex','FontSize',PlotSize,'FontName','Times');
set(gca,'FontSize',PlotSize,'Fontname','Times');
LgD=legend('Ionisation',...
    'Excitation',...
    'Heating');
set(LgD,'Interpreter','latex','Location','northwest','FontSize',PlotSize)
title('$m_{\rm{bh}} = 10^3 m_{\odot}$','Interpreter','latex','FontSize',PlotSize)
axis([10 2000 1E-2 1]);
