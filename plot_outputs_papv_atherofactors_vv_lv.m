% plot mass flux from vasa vasorum vs. into lymphatics

clear
close all

EXPORT = 1==1;
homepath = 'C:\Users\Willy\Work\PhD\ATLO\clean\';

%% load parameters
param = declare_fixed_parameters_1D(1==0);

load('output_1D_healthy_fluid.mat');
H.p_apv = p_apv;
for i = 1:numel(varpar)
    H.VV(i) = varpar(i).l_vv*10.*varpar(i).n_vv*param.l_ex*pi*param.d_vv; % cm/s/mmHg * mm^-2 * mm = /s/mmHg 
    H.VV(i) = H.VV(i)*(param.p_vv-param.sigma*(param.pi_vv-param.pi_ref)); % s^-1 (-> multiply by a reference pressure differential at zero hydrostatic tissue pressure)
    H.LV(i) = varpar(i).n_lv*1e6.*varpar(i).q_lv; % mm^-2 * m^2/s = s^-1
end

load('output_1D_athero_fluid.mat');
A.p_apv = p_apv;
for i = 1:numel(varpar)
    A.VV(i) = varpar(i).l_vv*10.*varpar(i).n_vv*param.l_ex*pi*param.d_vv; % cm/s/mmHg * mm^-2 * mm = /s/mmHg 
    A.VV(i) = A.VV(i)*(param.p_vv-param.sigma*(param.pi_vv-param.pi_ref)); % s^-1 (-> multiply by a reference pressure differential at zero hydrostatic tissue pressure)
    A.LV(i) = varpar(i).n_lv*1e6.*varpar(i).q_lv; % mm^-2 * m^2/s = s^-1
end


%% plot
figure('position',[50 50 800 800],'color','w');

labels = get_subplot_labels('a':'z',8);

ax1 = axes('Position',[0.32 0.56 0.4 0.4]);hold all;
set(gca,'fontsize',15,'TickLabelInterpreter','latex');
f = fill([0.6 1.4 1.4 0.6],[-10 -10 5 5],[0 0.8 0],'facealpha',0.5,'edgecolor','k','lines','--');
distributionPlot([H.p_apv',A.p_apv'],'xValues',[1 2],'color',0.75*[1 1 1],'distWidth',0.75,'showMM',6,'histOpt',1.1);
l_m = line([0 1],[50 50],'color','k','linew',2);
l_q = line([0 1],[50 50],'color','k','linew',1);
ylabel('$p_{\rm apv}$ (mmHg)','Interpreter','latex');
ylim([-20 20]);
xlim([0.6 3.5]); xticks([1 2]); xticklabels({'healthy';'athero.'})
legend([f l_m l_q],{'physiologic range';'median';'$1^{\rm st}/3^{\rm rd}$ quartile'},'location','southeast','fontsize',13,'Interpreter','latex');
text(-0.2,1,labels{1},'Units','normalized','FontSize',13,'fontweight','bold','FontName','times');

ax2 = axes('Position',[0.11 0.11 0.35 0.35]);hold all;
set(gca,'fontsize',15,'TickLabelInterpreter','latex');
ylim([0 2.5e-4])
ylabel('flux into lymphatics ($s^{-1}$)','Interpreter','latex')
distributionPlot([H.LV' A.LV'],'xvalues',[1 2],'color',0.75*[1 1 1],'distWidth',0.75,'showMM',6,'histOpt',1.1);
xlim([0.5 2.5])
xticklabels({'healthy';'athero.'})
text(-0.2,1.05,labels{2},'Units','normalized','FontSize',13,'fontweight','bold','FontName','times');

ax3 = axes('Position',[0.61 0.11 0.35 0.35]);hold all;
set(gca,'fontsize',15,'TickLabelInterpreter','latex');
ylim([0 2.5e-4])
ylabel({'flux from vasa vasorum';'at $p = 0\,(s^{-1})$'},'Interpreter','latex')
distributionPlot([H.VV' A.VV'],'xvalues',[1 2],'color',0.75*[1 1 1],'distWidth',0.75,'showMM',6,'histOpt',1.1);
xlim([0.5 2.5])
xticklabels({'healthy';'athero.'})
text(-0.2,1.05,labels{3},'Units','normalized','FontSize',13,'fontweight','bold','FontName','times');

if EXPORT
    exportgraphics(gcf,[homepath 'figures\output_p_apv_atherofactors_vv_lv.png'],'Resolution',300);
end
