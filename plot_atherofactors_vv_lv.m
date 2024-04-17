% plot mass flux from vasa vasorum vs. into lymphatics

clear
close all

EXPORT = 1==1;
homepath = 'C:\Users\Willy\Work\PhD\ATLO\clean\';

%% load parameters
param = declare_fixed_parameters_1D(1==0);

load('output_1D_healthy_fluid.mat');
for i = 1:numel(varpar)
    VV_H(i) = varpar(i).l_vv*10.*varpar(i).n_vv*param.l_ex*pi*param.d_vv; % cm/s/mmHg * mm^-2 * mm = /s/mmHg 
    VV_H(i) = VV_H(i)*(param.p_vv-param.sigma*(param.pi_vv-param.pi_ref)); % s^-1 (-> multiply by a reference pressure differential at zero hydrostatic tissue pressure)
    LV_H(i) = varpar(i).n_lv*1e6.*varpar(i).q_lv; % mm^-2 * m^2/s = s^-1
end

load('output_1D_athero_fluid.mat');
for i = 1:numel(varpar)
    VV_A(i) = varpar(i).l_vv*10.*varpar(i).n_vv*param.l_ex*pi*param.d_vv; % cm/s/mmHg * mm^-2 * mm = /s/mmHg 
    VV_A(i) = VV_A(i)*(param.p_vv-param.sigma*(param.pi_vv-param.pi_ref)); % s^-1 (-> multiply by a reference pressure differential at zero hydrostatic tissue pressure)
    LV_A(i) = varpar(i).n_lv*1e6.*varpar(i).q_lv; % mm^-2 * m^2/s = s^-1
end


%% plot
figure('position',[50 50 800 400],'color','w');
tiledlayout(1,2);
labels = get_subplot_labels('a':'z',8);

nexttile;hold all;
set(gca,'fontsize',13,'TickLabelInterpreter','latex');
ylim([0 2.5e-4])
ylabel('flux into lymphatics ($s^{-1}$)','Interpreter','latex')
distributionPlot([LV_H' LV_A'],'xvalues',[1 2],'color',0.75*[1 1 1],'distWidth',0.75,'showMM',6,'histOpt',1.1);
xlim([0.5 2.5])
xticklabels({'healthy';'atherosclerotic'})
text(-0.15,1.05,labels{1},'Units','normalized','FontSize',13,'fontweight','bold','FontName','times');

nexttile;hold all;
set(gca,'fontsize',13,'TickLabelInterpreter','latex');
ylim([0 2.5e-4])
ylabel('flux from vasa vasorum at $p = 0\,(s^{-1})$','Interpreter','latex')
distributionPlot([VV_H' VV_A'],'xvalues',[1 2],'color',0.75*[1 1 1],'distWidth',0.75,'showMM',6,'histOpt',1.1);
xlim([0.5 2.5])
xticklabels({'healthy';'atherosclerotic'})
text(-0.15,1.05,labels{2},'Units','normalized','FontSize',13,'fontweight','bold','FontName','times');

if EXPORT
    exportgraphics(gcf,[homepath 'figures\comparison_atherofactors_vv_lv.png'],'Resolution',300);
end
