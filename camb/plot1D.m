% An example for pretty 1D plot
clear
LineWidth=2
PlotSize=15
X_min=0
X_max=2*pi
Y_min=-1
Y_max=1

clc
% -- Get Array
x=0:0.01:2*pi;
y1=sin(x);
y2=cos(x);

% -- Get plot
clf
plot(x,y1,'k','LineWidth',LineWidth);hold on
plot(x,y2,'r','LineWidth',LineWidth);hold on
xlabel('$x$','Interpreter','latex','FontSize',PlotSize,'FontName','Times');
ylabel('$f(x)$','Interpreter','latex','FontSize',PlotSize,'FontName','Times');
set(gca,'FontSize',PlotSize,'Fontname','Times');
title('1D Example','Interpreter','latex','FontSize',PlotSize)
axis([X_min X_max Y_min Y_max]);
LgD=legend('$sin(x)$',...
    '$cos(x)$');
set(LgD,'Interpreter','latex','Location','best','FontSize',PlotSize)
% TICKS=10.^[-30:30];
% xticks([TICKS])
% yticks([TICKS])