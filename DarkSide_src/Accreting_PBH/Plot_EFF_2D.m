cd ~/cloud/GitHub/Dark_CosmoMC/DarkSide_src/Accreting_PBH/
clear

PlotSize=20
cid = 4

clc
load data/EFF.mat
f(:,:) = EFF(cid,:,:);
if cid == 1
    Title = 'Ionisation';
elseif cid == 3
    Title = 'Excitation';
elseif cid == 4
    Title = 'Heating';
end

clf
h=pcolor(M_Axis,zp_Axis,f);hold on
set(h, 'EdgeColor', 'none');
colormap('jet')
set(gca,'xdir','normal')
set(gca,'ydir','normal')
set(gca,'XScale','log') 
set(gca,'YScale','log') 
xlabel('$M_{\rm{bh}}\ [{\rm{M}}_{\odot}]$','Interpreter','latex','FontSize',PlotSize,'FontName','Times');
ylabel('$z$','Interpreter','latex','FontSize',PlotSize,'FontName','Times');
set(gca,'FontSize',PlotSize,'Fontname','Times');
xticks(10.^(-20:1:20))
CB=colorbar;
%set(get(CB,'Title'),'string','$L_{\rm{bh}}$','Interpreter','latex');
title(Title,'Interpreter','latex','FontSize',PlotSize)
