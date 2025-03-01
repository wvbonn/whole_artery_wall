clear
close all

currdir = pwd;
k = strfind(currdir,'\');
homepath = currdir(1:k(end));

EXPORT = 1==1;
labels = get_subplot_labels('a':'z',8);
addpath([homepath 'lumped\']);
col = [[0 0.6 0];[1 0 0.7]];

param_H = declare_fixed_parameters(1==0);
load("output_healthy_fluid.mat");
H.u_i = u_i; H.p_ia = p_ia; H.p_ap = p_ap; H.Q_v_l = Q_v_l; H.Dil = d_dil;
% get min and max PVAT thickness
for i = 1:numel(varpar); H.t_p(i) = varpar(i).t_p; end
H.r_ext = param_H.r_ia + param_H.t_a + [min(H.t_p) max(H.t_p)];

param_A = declare_fixed_parameters(1==1);
load('output_athero_fluid.mat');
A.u_i = u_i; A.p_ia = p_ia; A.p_ap = p_ap; A.Q_v_l = Q_v_l; A.Dil = d_dil;
% get min and max PVAT thickness
for i = 1:numel(varpar); A.t_p(i) = varpar(i).t_p; end
A.r_ext = param_A.r_ia + param_A.t_a + [min(A.t_p) max(A.t_p)];


return;

figure('position',[50 0 1000 900],'color','w');
tiledlayout(2,2);


%% inner-layer velocity

nexttile;hold all;
set(gca,'fontsize',15,'TickLabelInterpreter','latex');
f = fill([0.6 2.5 2.5 0.6],[0.025 0.025 0.055 0.055],[1 1 0],'facealpha',0.5,'edgecolor','k','lines','--');
distributionPlot([H.u_i',A.u_i'],'xValues',[1 2],'color',0.75*[1 1 1],'distWidth',0.75,'showMM',6,'histOpt',1.1);
l_m = line([0 1],[1 1],'color','k','linew',2);
l_q = line([0 1],[1 1],'color','k','linew',1);
ylabel('$\bar{u}_i (\mu{\rm m/s})$','Interpreter','latex'); 
ylim([0 0.09]); yticks(0:0.02:0.08);
xlim([0.6 3.8]); xticks([1 2]); xticklabels({'healthy';'athero.'});
legend([f l_m l_q],{'animal experiments';'median';'$1^{\rm st}/3^{\rm rd}$ quartile'},'location','northeast','fontsize',11,'Interpreter','latex');
text(-0.15,1,labels{1},'Units','normalized','FontSize',13,'fontweight','bold','FontName','times');


%% medial-adventitial pressure

nexttile;hold all;
set(gca,'fontsize',15,'TickLabelInterpreter','latex');
h = distributionPlot([H.p_ia',A.p_ia'],'xValues',[1 2],'color',0.75*[1 1 1],'distWidth',0.75,'showMM',6,'histOpt',1.1);
l_m = line([0 1],[50 50],'color','k','linew',2);
l_q = line([0 1],[50 50],'color','k','linew',1);
ylabel('$p_{ia}$ (mmHg)','Interpreter','latex'); 
ylim([-15 40]);
xlim([0.6 3.5]); xticks([1 2]); xticklabels({'healthy';'athero.'})
legend([l_m l_q],{'median';'$1^{\rm st}/3^{\rm rd}$ quartile'},'location','southeast','fontsize',11,'Interpreter','latex');
text(-0.15,1,labels{2},'Units','normalized','FontSize',13,'fontweight','bold','FontName','times');


% %% outer-adventitial pressure
% 
% nexttile;hold all;
% set(gca,'fontsize',15,'TickLabelInterpreter','latex');
% f = fill([0.6 1.4 1.4 0.6],[-10 -10 5 5],[0 0.8 0],'facealpha',0.5,'edgecolor','k','lines','--');
% distributionPlot([H.p_ap',A.p_ap'],'xValues',[1 2],'color',0.75*[1 1 1],'distWidth',0.75,'showMM',6,'histOpt',1.1);
% l_m = line([0 1],[50 50],'color','k','linew',2);
% l_q = line([0 1],[50 50],'color','k','linew',1);
% ylabel('$p_{\rm ap}$ (mmHg)','Interpreter','latex');
% ylim([-20 20]);
% xlim([0.6 3.5]); xticks([1 2]); xticklabels({'healthy';'athero.'})
% legend([f l_m l_q],{'physiologic range';'median';'$1^{\rm st}/3^{\rm rd}$ quartile'},'location','southeast','fontsize',13,'Interpreter','latex');
% text(-0.15,1,labels{1},'Units','normalized','FontSize',13,'fontweight','bold','FontName','times');
% 

%% share of vasa fluid in lymph

nexttile;hold all;
set(gca,'fontsize',15,'TickLabelInterpreter','latex');
ylabel('$Q_v/Q_{\ell}$','Interpreter','latex'); % yl = get(gca,'ylim'); % that line for conduits
l_m = line([0 1],[50 50],'color','k','linew',2);
l_q = line([0 1],[50 50],'color','k','linew',1);
ylim([0 1]);
yticks([0 0.2 0.4 0.6 0.8 1]);
xlim([0.6 3.5]);
distributionPlot([H.Q_v_l',A.Q_v_l'],'xValues',[1 2],'color',0.75*[1 1 1],'distWidth',0.75,'showMM',6,'histOpt',1.1);
xlim([0.6 3.8]); xticks([1 2]); xticklabels({'healthy';'athero.'})
legend([l_m l_q],{'median';'$1^{\rm st}/3^{\rm rd}$ quartile'},'location','southeast','fontsize',11,'Interpreter','latex');
text(-0.15,1,labels{3},'Units','normalized','FontSize',13,'fontweight','bold','FontName','times');



%% dilution distances
r_ia = param_H.r_ia;

H.r_lu = r_ia-param_H.t_i;
H.r_ap = r_ia+param_H.t_a;

A.r_lu = r_ia-param_A.t_i;
A.r_ap = r_ia+param_A.t_a;


nexttile;
hold all;
set(gca,'fontsize',15,'TickLabelInterpreter','latex');
f_inH = fill([0.6 1.4 1.4 0.6],[H.r_lu*[1 1] r_ia*[1 1]],[0 0 1],'linew',1.5,'facealpha',1,'edgecolor','none');
f_inA = fill([1.6 2.4 2.4 1.6],[A.r_lu*[1 1] r_ia*[1 1]],[0 0 1],'linew',1.5,'facealpha',1,'edgecolor','none');
f_eelH = fill([0.6 1.4 1.4 0.6],[r_ia*[1 1] H.r_ap*[1 1]],[0.1 0 1],'linew',1.5,'facealpha',0.4,'edgecolor','none');
f_eelA = fill([1.6 2.4 2.4 1.6],[r_ia*[1 1] A.r_ap*[1 1]],[0.1 0 1],'linew',1.5,'facealpha',0.4,'edgecolor','none');
f_pvHmin = fill([0.6 1.4 1.4 0.6],[H.r_ap*[1 1] H.r_ext(1)*[1 1]],[0.4 0 1],'linew',1.5,'facealpha',0.2,'edgecolor','none');
f_pvAmin = fill([1.6 2.4 2.4 1.6],[A.r_ap*[1 1] A.r_ext(1)*[1 1]],[0.4 0 1],'linew',1.5,'facealpha',0.2,'edgecolor','none');
f_pvHmax = fill([0.6 1.4 1.4 0.6],[H.r_ext(1)*[1 1] H.r_ext(2)*[1 1]],[1 1 1],'linew',1.5,'lines','--','edgecolor',[1 0.75 1]);
f_pvAmax = fill([1.6 2.4 2.4 1.6],[A.r_ext(1)*[1 1] A.r_ext(2)*[1 1]],[1 1 1],'linew',1.5,'lines','--','edgecolor',[1 0.75 1]);
distributionPlot([H.Dil',A.Dil'],'xValues',[1 2],'color',0.75*[1 1 1],'distWidth',0.75,'showMM',6,'histOpt',1.1);
ylabel('dilution distance (mm)','Interpreter','latex'); 
ylim([1 4.68])
xlim([0.6 4]); xticks([1 2]); xticklabels({'healthy';'athero.'});
legend([f_pvAmax f_pvAmin f_eelA f_inA],'max. PVAT','min. PVAT','adventitia','inner layers','location','southeast','fontsize',11,'Interpreter','latex');
text(-0.15,1,labels{4},'Units','normalized','FontSize',13,'fontweight','bold','FontName','times');

if EXPORT
    exportgraphics(gcf,[homepath 'figures\output_ui_pia_qvl_ddil.png'],'Resolution',300);
end



