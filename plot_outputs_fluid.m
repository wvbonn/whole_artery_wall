clear
close all

EXPORT = 1==0;
GSA = 1==0;
homepath = 'C:\Users\Willy\Work\PhD\ATLO\clean\';
labels = get_subplot_labels('a':'z',8);
addpath([homepath 'lumped\']);
col = [[0 0.6 0];[1 0 0.7]];

if GSA
    load([homepath 'lumped\output_GSA_0D_healthy.mat']);
    Sglob(Sglob<0) = 0;
    Hv.t = Sglob(1,1:8); %Hv = sortmaincoeffs(Hv);
    Hpe.t = Sglob(2,1:8);
    Hpa.t = Sglob(3,1:8); %Hp = sortmaincoeffs(Hp);
    Hqv.t = Sglob(4,1:8); %Hqv = sortmaincoeffs(Hqv);
    
    load([homepath 'lumped\output_GSA_0D_athero.mat']);
    Sglob(Sglob<0) = 0;
    Av.t = Sglob(1,1:8); %Hv = sortmaincoeffs(Hv); [1:8 10:12]
    Ape.t = Sglob(2,1:8);
    Apa.t = Sglob(3,1:8); %Hp = sortmaincoeffs(Hp);
    Aqv.t = Sglob(4,1:8); %Hqv = sortmaincoeffs(Hqv);
end

xticknames_H = {'$k_{\rm im}$';'$k_{\rm a}$';'$t_{\rm pv}$';'$p_{\rm lum}$';'$l_{\rm p,vv}$';'$N_{\rm vv}$';'$N_{\rm lv}$';'$\tilde{q}_{lv}$'};
xticknames_A = xticknames_H;

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

%% Dilution distances
r_eel = param_H.r_eel;

H.r_lu = r_eel-param_H.t_im;
H.r_apv = r_eel+param_H.t_a;

A.r_lu = r_eel-param_A.t_im;
A.r_apv = r_eel+param_A.t_a;

figure('position',[50 100 500 400],'color','w');
% tiledlayout(1,2,'TileSpacing','compact');
% labels = get_subplot_labels('A':'Z',6);
% 
% nexttile;
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
xlim([0.6 4]); xticks([1 2]); xticklabels({'healthy';'atherosclerotic'});
legend([f_pvAmax f_pvAmin f_eelA f_inA],'max. PVAT','min. PVAT','adventitia','inner layers','location','southeast','fontsize',13,'Interpreter','latex');


if EXPORT
    exportgraphics(gcf,[homepath 'figures\dilution_distances.png'],'Resolution',300);
end


%% inner-layer velocity


figure('Position',[50 50 1000 450],'color','w');
tiledlayout(2,2,'TileSpacing','compact');

nexttile([2 1]);hold all;
set(gca,'fontsize',15,'TickLabelInterpreter','latex');
f = fill([0.6 2.5 2.5 0.6],[0.025 0.025 0.055 0.055],[1 1 0],'facealpha',0.5,'edgecolor','k','lines','--');
distributionPlot([H.v_im',A.v_im'],'xValues',[1 2],'color',0.75*[1 1 1],'distWidth',0.75,'showMM',6,'histOpt',1.1);
l_m = line([0 1],[1 1],'color','k','linew',2);
l_q = line([0 1],[1 1],'color','k','linew',1);
ylabel('$\bar{v}_{\rm im} (\mu{\rm m/s})$','Interpreter','latex'); 
ylim([0 0.09]); yticks(0:0.02:0.08);
xlim([0.6 3.8]); xticks([1 2]); xticklabels({'healthy';'atherosclerotic'});
legend([f l_m l_q],{'animal experiments';'median';'$1^{\rm st}/3^{\rm rd}$ quartile'},'location','northeast','fontsize',13,'Interpreter','latex');
text(-0.15,1,labels{1},'Units','normalized','FontSize',13,'fontweight','bold','FontName','times');

if GSA
    nexttile;
    bar(Hv.t,'facecolor',0.4*[1 1 1]);
    xticklabels(xticknames_H);
    xtickangle(45);
    hold on;box off;
    title('healthy','fontsize',14,'Interpreter','latex');
    set(gca,'fontsize',14,'TickLabelInterpreter','latex');
    ylim([0 0.9]);
    yticks([0.05 0.2 0.4 0.6 0.8]);
    ylabel('global sensitivity','Interpreter','latex');
    xl = get(gca,'XLim');
    line(xl,[0.05 0.05],'lines','--','color','k','linew',1.25);
    
    t = text(1,Hv.t(1)+0.05,'+','FontSize',15,'FontWeight','bold','FontName','times');
    t.Position(1) = 1-t.Extent(3)/2;
    t = text(2,Hv.t(2)+0.05,'+','FontSize',15,'FontWeight','bold','FontName','times');
    t.Position(1) = 2-t.Extent(3)/2;
    t = text(4,Hv.t(4)+0.05,'+','FontSize',15,'FontWeight','bold','FontName','times');
    t.Position(1) = 4-t.Extent(3)/2;
    text(-0.1,1,labels{2},'Units','normalized','FontSize',13,'fontweight','bold','FontName','times');
    
    nexttile;
    bar(Av.t,'facecolor',0.4*[1 1 1]);
    xticklabels(xticknames_A);
    xtickangle(45);
    hold on;box off;
    title('atherosclerotic','fontsize',14,'Interpreter','latex');
    set(gca,'fontsize',14,'TickLabelInterpreter','latex');
    ylim([0 0.9]);
    yticks([0.05 0.2 0.4 0.6 0.8]);
    ylabel('global sensitivity','Interpreter','latex');
    xl = get(gca,'XLim');
    line(xl,[0.05 0.05],'lines','--','color','k','linew',1.25);
    
    t = text(1,Av.t(1)+0.05,'+','FontSize',15,'FontWeight','bold','FontName','times');
    t.Position(1) = 1-t.Extent(3)/2;
    t = text(2,Av.t(2)+0.05,'+','FontSize',15,'FontWeight','bold','FontName','times');
    t.Position(1) = 2-t.Extent(3)/2;
    t = text(4,Av.t(4)+0.05,'+','FontSize',15,'FontWeight','bold','FontName','times');
    t.Position(1) = 4-t.Extent(3)/2;
    text(-0.1,1,labels{3},'Units','normalized','FontSize',13,'fontweight','bold','FontName','times');
end

if EXPORT
    exportgraphics(gcf,[homepath 'figures\output_v_im.png'],'Resolution',300);
end

%% medial-adventitial pressure


figure('Position',[50 50 1000 450],'color','w');
tiledlayout(2,2,'TileSpacing','compact');

nexttile([2 1]);hold all;
set(gca,'fontsize',15,'TickLabelInterpreter','latex');
h = distributionPlot([H.p_ma',A.p_ma'],'xValues',[1 2],'color',0.75*[1 1 1],'distWidth',0.75,'showMM',6,'histOpt',1.1);
l_m = line([0 1],[50 50],'color','k','linew',2);
l_q = line([0 1],[50 50],'color','k','linew',1);
ylabel('$p_{\rm ma}$ (mmHg)','Interpreter','latex'); 
ylim([-15 40]);
xlim([0.6 3.5]); xticks([1 2]); xticklabels({'healthy';'atherosclerotic'})
legend([l_m l_q],{'median';'$1^{\rm st}/3^{\rm rd}$ quartile'},'location','southeast','fontsize',13,'Interpreter','latex');
text(-0.15,1,labels{1},'Units','normalized','FontSize',13,'fontweight','bold','FontName','times');

if GSA
    nexttile;
    bar(Hpe.t,'facecolor',0.4*[1 1 1]);
    xticklabels(xticknames_H);
    xtickangle(45);
    hold on;box off;
    title('healthy','fontsize',14,'Interpreter','latex');
    set(gca,'fontsize',14,'TickLabelInterpreter','latex');
    ylim([0 0.9]);
    yticks([0.05 0.2 0.4 0.6 0.8]);
    ylabel('global sensitivity','Interpreter','latex');
    xl = get(gca,'XLim');
    line(xl,[0.05 0.05],'lines','--','color','k','linew',1.25);
    
    t = text(1,Hpe.t(1)+0.05,'+','FontSize',15,'FontWeight','bold','FontName','times');
    t.Position(1) = 1-t.Extent(3)/2;
    t = text(2,Hpe.t(2)+0.05,'-','FontSize',15,'FontWeight','bold','FontName','times');
    t.Position(1) = 2-t.Extent(3)/2;
    t = text(8,Hpe.t(8)+0.05,'-','FontSize',15,'FontWeight','bold','FontName','times');
    t.Position(1) = 8-t.Extent(3)/2;
    text(-0.1,1,labels{2},'Units','normalized','FontSize',13,'fontweight','bold','FontName','times');
    
    nexttile;
    bar(Ape.t,'facecolor',0.4*[1 1 1]);
    xticklabels(xticknames_A);
    xtickangle(45);
    hold on;box off;
    title('atherosclerotic','fontsize',14,'Interpreter','latex');
    set(gca,'fontsize',14,'TickLabelInterpreter','latex');
    ylim([0 0.9]);
    yticks([0.05 0.2 0.4 0.6 0.8]);
    ylabel('global sensitivity','Interpreter','latex');
    xl = get(gca,'XLim');
    line(xl,[0.05 0.05],'lines','--','color','k','linew',1.25);
    
    t = text(1,Ape.t(1)+0.05,'+','FontSize',15,'FontWeight','bold','FontName','times');
    t.Position(1) = 1-t.Extent(3)/2;
    t = text(2,Ape.t(2)+0.05,'-','FontSize',15,'FontWeight','bold','FontName','times');
    t.Position(1) = 2-t.Extent(3)/2;
    t = text(8,Ape.t(8)+0.05,'-','FontSize',15,'FontWeight','bold','FontName','times');
    t.Position(1) = 8-t.Extent(3)/2;
    text(-0.1,1,labels{3},'Units','normalized','FontSize',13,'fontweight','bold','FontName','times');
end

if EXPORT
    exportgraphics(gcf,[homepath 'figures\output_p_ma.png'],'Resolution',300);
end

%% outer-adventitial pressure

figure('Position',[50 50 1000 450],'color','w');
tiledlayout(2,2,'TileSpacing','compact');

nexttile([2 1]);hold all;
set(gca,'fontsize',15,'TickLabelInterpreter','latex');
f = fill([0.6 1.4 1.4 0.6],[-10 -10 5 5],[0 0.8 0],'facealpha',0.5,'edgecolor','k','lines','--');
distributionPlot([H.p_apv',A.p_apv'],'xValues',[1 2],'color',0.75*[1 1 1],'distWidth',0.75,'showMM',6,'histOpt',1.1);
l_m = line([0 1],[50 50],'color','k','linew',2);
l_q = line([0 1],[50 50],'color','k','linew',1);
ylabel('$p_{\rm apv}$ (mmHg)','Interpreter','latex');
ylim([-20 20]);
xlim([0.6 3.5]); xticks([1 2]); xticklabels({'healthy';'atherosclerotic'})
legend([f l_m l_q],{'physiologic range';'median';'$1^{\rm st}/3^{\rm rd}$ quartile'},'location','southeast','fontsize',13,'Interpreter','latex');
text(-0.15,1,labels{1},'Units','normalized','FontSize',13,'fontweight','bold','FontName','times');

if GSA
    nexttile;
    bar(Hpa.t,'facecolor',0.4*[1 1 1]);
    xticklabels(xticknames_H);
    xtickangle(45); 
    hold on;box off;
    title('healthy','fontsize',14,'Interpreter','latex');
    set(gca,'fontsize',14,'TickLabelInterpreter','latex');
    ylim([0 0.9]);
    yticks([0.05 0.2 0.4 0.6 0.8]);
    ylabel('global sensitivity','Interpreter','latex');
    xl = get(gca,'XLim');
    line(xl,[0.05 0.05],'lines','--','color','k','linew',1.25);    
    t = text(1,Hpa.t(1)+0.05,'+','FontSize',15,'FontWeight','bold','FontName','times');
    t.Position(1) = 1-t.Extent(3)/2;
    t = text(2,Hpa.t(2)+0.05,'+','FontSize',15,'FontWeight','bold','FontName','times');
    t.Position(1) = 2-t.Extent(3)/2;
    t = text(3,Hpa.t(3)+0.05,'-','FontSize',15,'FontWeight','bold','FontName','times');
    t.Position(1) = 3-t.Extent(3)/2;
    t = text(6,Hpa.t(6)+0.05,'+','FontSize',15,'FontWeight','bold','FontName','times');
    t.Position(1) = 6-t.Extent(3)/2;
    t = text(8,Hpa.t(8)+0.05,'-','FontSize',15,'FontWeight','bold','FontName','times');
    t.Position(1) = 8-t.Extent(3)/2;
    text(-0.1,1,labels{2},'Units','normalized','FontSize',13,'fontweight','bold','FontName','times');
    
    nexttile;
    bar(Apa.t,'facecolor',0.4*[1 1 1]);
    xticklabels(xticknames_A);
    xtickangle(45);
    hold on;box off;
    title('atherosclerotic','fontsize',14,'Interpreter','latex');
    set(gca,'fontsize',14,'TickLabelInterpreter','latex');
    ylim([0 0.9]);
    yticks([0.05 0.2 0.4 0.6 0.8]);
    ylabel('global sensitivity','Interpreter','latex');
    xl = get(gca,'XLim');
    line(xl,[0.05 0.05],'lines','--','color','k','linew',1.25);    
    t = text(1,Apa.t(1)+0.05,'+','FontSize',15,'FontWeight','bold','FontName','times');
    t.Position(1) = 1-t.Extent(3)/2;
    t = text(2,Apa.t(2)+0.05,'+','FontSize',15,'FontWeight','bold','FontName','times');
    t.Position(1) = 2-t.Extent(3)/2;
    t = text(3,Apa.t(3)+0.05,'-','FontSize',15,'FontWeight','bold','FontName','times');
    t.Position(1) = 3-t.Extent(3)/2;
    t = text(6,Apa.t(6)+0.05,'+','FontSize',15,'FontWeight','bold','FontName','times');
    t.Position(1) = 6-t.Extent(3)/2;
    t = text(8,Apa.t(8)+0.05,'-','FontSize',15,'FontWeight','bold','FontName','times');
    t.Position(1) = 8-t.Extent(3)/2;
    text(-0.1,1,labels{3},'Units','normalized','FontSize',11,'fontweight','bold','FontName','times');
end

if EXPORT
    exportgraphics(gcf,[homepath 'figures\output_p_apv.png'],'Resolution',300);
end


%% share of vasa fluid in lymph

figure('Position',[50 50 1000 450],'color','w');
tiledlayout(2,2,'TileSpacing','compact');

nexttile([2 1]);hold all;
set(gca,'fontsize',15,'TickLabelInterpreter','latex');
ylabel('$Q_{\rm vv}/Q_{\rm lv}$','Interpreter','latex'); % yl = get(gca,'ylim'); % that line for conduits
l_m = line([0 1],[50 50],'color','k','linew',2);
l_q = line([0 1],[50 50],'color','k','linew',1);
ylim([0 1]);
yticks([0 0.2 0.4 0.6 0.8 1]);
xlim([0.6 3.5]);
distributionPlot([H.Q_vv_lv',A.Q_vv_lv'],'xValues',[1 2],'color',0.75*[1 1 1],'distWidth',0.75,'showMM',6,'histOpt',1.1);
xlim([0.6 3.8]); xticks([1 2]); xticklabels({'healthy';'atherosclerotic'})
legend([l_m l_q],{'median';'$1^{\rm st}/3^{\rm rd}$ quartile'},'location','southeast','fontsize',13,'Interpreter','latex');
text(-0.15,1,labels{1},'Units','normalized','FontSize',13,'fontweight','bold','FontName','times');

if GSA
    nexttile;
    bar(Hqv.t,'facecolor',0.4*[1 1 1]);
    xticklabels(xticknames_H);
    xtickangle(45);
    hold on;box off;
    title('healthy','fontsize',14,'Interpreter','latex');
    set(gca,'fontsize',14,'TickLabelInterpreter','latex');
    ylim([0 0.9]);
    ylabel('global sensitivity','Interpreter','latex');
    yticks([0.05 0.2 0.4 0.6 0.8]);
    xl = get(gca,'XLim');
    line(xl,[0.05 0.05],'lines','--','color','k','linew',1.25);    
    t = text(1,Hqv.t(1)+0.05,'-','FontSize',15,'FontWeight','bold','FontName','times');
    t.Position(1) = 1-t.Extent(3)/2;
    t = text(2,Hqv.t(2)+0.05,'-','FontSize',15,'FontWeight','bold','FontName','times');
    t.Position(1) = 2-t.Extent(3)/2;
    t = text(5,Hqv.t(5)+0.05,'+','FontSize',15,'FontWeight','bold','FontName','times');
    t.Position(1) = 5-t.Extent(3)/2;
    t = text(6,Hqv.t(6)+0.05,'+','FontSize',15,'FontWeight','bold','FontName','times');
    t.Position(1) = 6-t.Extent(3)/2;
    t = text(8,Hqv.t(8)+0.05,'+','FontSize',15,'FontWeight','bold','FontName','times');
    t.Position(1) = 8-t.Extent(3)/2;
    text(-0.1,1,labels{2},'Units','normalized','FontSize',13,'fontweight','bold','FontName','times');
    
    nexttile;
    bar(Aqv.t,'facecolor',0.4*[1 1 1]);
    xticklabels(xticknames_A);
    xtickangle(45);
    hold on;box off;
    title('atherosclerotic','fontsize',14,'Interpreter','latex');
    set(gca,'fontsize',14,'TickLabelInterpreter','latex');
    ylim([0 0.9]);
    ylabel('global sensitivity','Interpreter','latex');
    yticks([0.05 0.2 0.4 0.6 0.8]);
    xl = get(gca,'XLim');
    line(xl,[0.05 0.05],'lines','--','color','k','linew',1.25);    
    t = text(1,Aqv.t(1)+0.05,'-','FontSize',15,'FontWeight','bold','FontName','times');
    t.Position(1) = 1-t.Extent(3)/2;
    t = text(2,Aqv.t(2)+0.05,'-','FontSize',15,'FontWeight','bold','FontName','times');
    t.Position(1) = 2-t.Extent(3)/2;
    t = text(8,Aqv.t(8)+0.05,'+','FontSize',15,'FontWeight','bold','FontName','times');
    t.Position(1) = 8-t.Extent(3)/2;
    text(-0.1,1,labels{3},'Units','normalized','FontSize',13,'fontweight','bold','FontName','times');
end

if EXPORT
    exportgraphics(gcf,[homepath 'figures\output_Q_vv_lv.png'],'Resolution',300);
end



