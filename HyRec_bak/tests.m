cd /Users/cangtao/IHEPBox/Projects/GitHub/Dark_CosmoMC/HyRec
clear

LineWidth=2
PlotSize=20

clc

d=load('tmp.dat');
zp=d(:,1);
r=d(:,2);

! cp ~/FUNCTIONS/CurveReader.m ./
File = '/Users/cangtao/OneDrive/CMB/Result/CrossCheck/0709.0524.Fig3.Left.Top.txt';
File = '/Users/cangtao/OneDrive/CMB/Result/CrossCheck/0709.0524.Fig3.Right.Top.txt';
[x,y]=CurveReader(File,100,1);
! rm CurveReader.m
zp2=10.^x;
r2=10.^y;

clf
loglog(zp,r,'k','LineWidth',LineWidth);hold on
semilogx(zp2,r2,'r','LineWidth',LineWidth);hold on
xlabel('$z$','Interpreter','latex','FontSize',PlotSize,'FontName','Times');
ylabel('$\dot{m}$','Interpreter','latex','FontSize',PlotSize,'FontName','Times');
set(gca,'FontSize',PlotSize,'Fontname','Times');
LgD=legend('Mine',...
    'Paper');
set(LgD,'Interpreter','latex','Location','best','FontSize',PlotSize)
