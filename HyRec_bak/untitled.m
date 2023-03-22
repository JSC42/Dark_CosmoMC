cd /Users/cangtao/IHEPBox/Projects/GitHub/Dark_CosmoMC/HyRec
clear

LineWidth=2
PlotSize=20

! echo 1E3 > tmp_in.dat
% F = '/Users/cangtao/OneDrive/CMB/Result/CrossCheck/2303.06616.fig1.red_solid.txt'
% F = '/Users/cangtao/OneDrive/CMB/Result/CrossCheck/2303.06616.fig1.green_dashed.txt'
F = '/Users/cangtao/OneDrive/CMB/Result/CrossCheck/2303.06616.fig1.pink_dashed.txt'
% F = '/Users/cangtao/OneDrive/CMB/Result/CrossCheck/2303.06616.fig1.blue_dashed.txt'

clc
! cp ~/FUNCTIONS/CurveReader.m ./
[x,y]=CurveReader(F,100,1);
! rm CurveReader.m
zp1 = 10.^x+1;
md1 = 10.^y;

%d=load('tmp.dat');
! ./a.out < tmp_in.dat > tmp_out.dat
d=load('tmp_out.dat');
zp2=d(:,1) + 1;
md2=d(:,2);

clf
semilogx(zp1,md1,'k','LineWidth',LineWidth);hold on
semilogx(zp2,md2,'r','LineWidth',LineWidth);hold on
xlabel('$z$','Interpreter','latex','FontSize',PlotSize,'FontName','Times');
ylabel('$y$','Interpreter','latex','FontSize',PlotSize,'FontName','Times');
set(gca,'FontSize',PlotSize,'Fontname','Times');
%axis([-inf inf -inf inf])
