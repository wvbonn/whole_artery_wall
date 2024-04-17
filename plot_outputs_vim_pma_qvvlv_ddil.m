clear
close all

EXPORT = 1==1;
homepath = 'C:\Users\Willy\Work\PhD\ATLO\clean\';
labels = get_subplot_labels('a':'z',8);
addpath([homepath 'lumped\']);
col = [[0 0.6 0];[1 0 0.7]];

param_H = declare_fixed_parameters_1D(1==0);
load("output_1D_healthy_fluid.mat");
H.v_im = v_im; H.p_ma = p_ma; H.p_apv = p_apv; H.Q_vv_lv = Q_vv_lv; H.Dil = d_dil;
% get min and max PVAT thickness
for i = 1:numel(varpar); H.t_pv(i) = varpar(i).t_pv; end
H.r_pvext = param_H.r_eel + param_H.t_a + [min(H.t_pv) max(H.t_pv)];

param_A = declare_fixed_parameters_1D(1==1);
load('output_1D_athero_fluid.mat');
A.v_im = v_im; A.p_ma = p_ma; A.p_apv = p_apv; A.Q_vv_lv = Q_vv_lv; A.Dil = d_dil;
% get min and max PVAT thickness
for i = 1:numel(varpar); A.t_pv(i) = varpar(i).t_pv; end
A.r_pvext = param_A.r_eel + param_A.t_a + [min(A.t_pv) max(A.t_pv)];


figure('position',[50 0 1000 900],'color','w');
tiledlayout(2,2);


%% inner-layer velocity

nexttile;hold all;
set(gca,'fontsize',15,'TickLabelInterpreter','latex');
f = fill([0.6 2.5 2.5 0.6],[0.025 0.025 0.055 0.055],[1 1 0],'facealpha',0.5,'edgecolor','k','lines','--');
distributionPlot([H.v_im',A.v_im'],'xValues',[1 2],'color',0.75*[1 1 1],'distWidth',0.75,'showMM',6,'histOpt',1.1);
l_m = line([0 1],[1 1],'color','k','linew',2);
l_q = line([0 1],[1 1],'color','k','linew',1);
ylabel('$\bar{v}_{\rm im} (\mu{\rm m/s})$','Interpreter','latex'); 
ylim([0 0.09]); yticks(0:0.02:0.08);
xlim([0.6 3.8]); xticks([1 2]); xticklabels({'healthy';'athero.'});
legend([f l_m l_q],{'animal experiments';'median';'$1^{\rm st}/3^{\rm rd}$ quartile'},'location','northeast','fontsize',11,'Interpreter','latex');
text(-0.15,1,labels{1},'Units','normalized','FontSize',13,'fontweight','bold','FontName','times');


%% medial-adventitial pressure

nexttile;hold all;
set(gca,'fontsize',15,'TickLabelInterpreter','latex');
h = distributionPlot([H.p_ma',A.p_ma'],'xValues',[1 2],'color',0.75*[1 1 1],'distWidth',0.75,'showMM',6,'histOpt',1.1);
l_m = line([0 1],[50 50],'color','k','linew',2);
l_q = line([0 1],[50 50],'color','k','linew',1);
ylabel('$p_{\rm ma}$ (mmHg)','Interpreter','latex'); 
ylim([-15 40]);
xlim([0.6 3.5]); xticks([1 2]); xticklabels({'healthy';'athero.'})
legend([l_m l_q],{'median';'$1^{\rm st}/3^{\rm rd}$ quartile'},'location','southeast','fontsize',11,'Interpreter','latex');
text(-0.15,1,labels{2},'Units','normalized','FontSize',13,'fontweight','bold','FontName','times');


% %% outer-adventitial pressure
% 
% nexttile;hold all;
% set(gca,'fontsize',15,'TickLabelInterpreter','latex');
% f = fill([0.6 1.4 1.4 0.6],[-10 -10 5 5],[0 0.8 0],'facealpha',0.5,'edgecolor','k','lines','--');
% distributionPlot([H.p_apv',A.p_apv'],'xValues',[1 2],'color',0.75*[1 1 1],'distWidth',0.75,'showMM',6,'histOpt',1.1);
% l_m = line([0 1],[50 50],'color','k','linew',2);
% l_q = line([0 1],[50 50],'color','k','linew',1);
% ylabel('$p_{\rm apv}$ (mmHg)','Interpreter','latex');
% ylim([-20 20]);
% xlim([0.6 3.5]); xticks([1 2]); xticklabels({'healthy';'athero.'})
% legend([f l_m l_q],{'physiologic range';'median';'$1^{\rm st}/3^{\rm rd}$ quartile'},'location','southeast','fontsize',13,'Interpreter','latex');
% text(-0.15,1,labels{1},'Units','normalized','FontSize',13,'fontweight','bold','FontName','times');
% 

%% share of vasa fluid in lymph

nexttile;hold all;
set(gca,'fontsize',15,'TickLabelInterpreter','latex');
ylabel('$Q_{\rm vv}/Q_{\rm lv}$','Interpreter','latex'); % yl = get(gca,'ylim'); % that line for conduits
l_m = line([0 1],[50 50],'color','k','linew',2);
l_q = line([0 1],[50 50],'color','k','linew',1);
ylim([0 1]);
yticks([0 0.2 0.4 0.6 0.8 1]);
xlim([0.6 3.5]);
distributionPlot([H.Q_vv_lv',A.Q_vv_lv'],'xValues',[1 2],'color',0.75*[1 1 1],'distWidth',0.75,'showMM',6,'histOpt',1.1);
xlim([0.6 3.8]); xticks([1 2]); xticklabels({'healthy';'athero.'})
legend([l_m l_q],{'median';'$1^{\rm st}/3^{\rm rd}$ quartile'},'location','southeast','fontsize',11,'Interpreter','latex');
text(-0.15,1,labels{3},'Units','normalized','FontSize',13,'fontweight','bold','FontName','times');



%% dilution distances
r_eel = param_H.r_eel;

H.r_lu = r_eel-param_H.t_im;
H.r_apv = r_eel+param_H.t_a;

A.r_lu = r_eel-param_A.t_im;
A.r_apv = r_eel+param_A.t_a;


nexttile;
hold all;
set(gca,'fontsize',15,'TickLabelInterpreter','latex');
f_inH = fill([0.6 1.4 1.4 0.6],[H.r_lu*[1 1] r_eel*[1 1]],[0 0 1],'linew',1.5,'facealpha',1,'edgecolor','none');
f_inA = fill([1.6 2.4 2.4 1.6],[A.r_lu*[1 1] r_eel*[1 1]],[0 0 1],'linew',1.5,'facealpha',1,'edgecolor','none');
f_eelH = fill([0.6 1.4 1.4 0.6],[r_eel*[1 1] H.r_apv*[1 1]],[0.1 0 1],'linew',1.5,'facealpha',0.4,'edgecolor','none');
f_eelA = fill([1.6 2.4 2.4 1.6],[r_eel*[1 1] A.r_apv*[1 1]],[0.1 0 1],'linew',1.5,'facealpha',0.4,'edgecolor','none');
f_pvHmin = fill([0.6 1.4 1.4 0.6],[H.r_apv*[1 1] H.r_pvext(1)*[1 1]],[0.4 0 1],'linew',1.5,'facealpha',0.2,'edgecolor','none');
f_pvAmin = fill([1.6 2.4 2.4 1.6],[A.r_apv*[1 1] A.r_pvext(1)*[1 1]],[0.4 0 1],'linew',1.5,'facealpha',0.2,'edgecolor','none');
f_pvHmax = fill([0.6 1.4 1.4 0.6],[H.r_pvext(1)*[1 1] H.r_pvext(2)*[1 1]],[1 1 1],'linew',1.5,'lines','--','edgecolor',[1 0.75 1]);
f_pvAmax = fill([1.6 2.4 2.4 1.6],[A.r_pvext(1)*[1 1] A.r_pvext(2)*[1 1]],[1 1 1],'linew',1.5,'lines','--','edgecolor',[1 0.75 1]);
distributionPlot([H.Dil',A.Dil'],'xValues',[1 2],'color',0.75*[1 1 1],'distWidth',0.75,'showMM',6,'histOpt',1.1);
ylabel('dilution distance (mm)','Interpreter','latex'); 
ylim([1 4.5])
xlim([0.6 4]); xticks([1 2]); xticklabels({'healthy';'athero.'});
legend([f_pvAmax f_pvAmin f_eelA f_inA],'max. PVAT','min. PVAT','adventitia','inner layers','location','southeast','fontsize',11,'Interpreter','latex');
text(-0.15,1,labels{4},'Units','normalized','FontSize',13,'fontweight','bold','FontName','times');

if EXPORT
    exportgraphics(gcf,[homepath 'figures\output_vim_pma_qvvlv_ddil.png'],'Resolution',300);
end



