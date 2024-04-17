% comparison of p_apv between 1D model and 0D model. Direct comparison
% between cases with same parameters

clear
close all

EXPORT = 1==1;
homepath = 'C:\Users\Willy\Work\PhD\ATLO\clean\';

figure('position',[50 50 800 400],'color','w');
tiledlayout(1,2,'TileSpacing','compact');
labels = get_subplot_labels('a':'z',6);

%% healthy
ATHERO = 1==0;
load([homepath 'radial\output_1D_healthy_fluid.mat']);

param = declare_fixed_parameters_0D(ATHERO);

for i = 1:size(p_apv,2)
    Out = lumped_model(varpar(i),param);
    p_0D(i) = Out(3);
    p_1D(i) = p_apv(i);
end

ax1 = nexttile;hold all;
set(gca,'FontSize',14,'TickLabelInterpreter','latex');
title('healthy','Interpreter','latex')
text(-0.2,1,labels{1},'Units','normalized','FontSize',13,'fontweight','bold','FontName','times');
xlabel('$p_{\rm apv}$ (radial)','Interpreter','latex')
ylabel('$p_{\rm apv}$ (lumped-parameter)','Interpreter','latex')
[p_1Dsort,sortIdx] = sort(p_1D);
pdiffh = p_1Dsort-p_0D(sortIdx);
f5 = fill([-20 -15 15 15 10 -20],[-20 -20 10 15 15 -15],[0 0.8 0],'FaceAlpha',0.2,'EdgeColor','none');
f2 = fill([-20 -18 15 15 13 -20],[-20 -20 13 15 15 -18],[0 0.8 0],'FaceAlpha',0.45,'EdgeColor','none');
plot([-20 15],[-20 15],'linew',1.5,'color','k');
scatter(p_1Dsort,p_0D(sortIdx),10,'filled','+','MarkerEdgeColor',[0 0 0.8]);
text(-15,12,'{\it paired comparisons}','FontSize',13,'Interpreter','latex');
xlim([-20 15])
ylim([-20 15])
legend([f2 f5],{'$\pm$ 2 mmHg';'$\pm$ 5 mmHg'},'fontsize',14,'Interpreter','latex','location','southeast');


%% athero
clear p_1D p_0D p_1Dsort

ATHERO = 1==1;
load([homepath 'radial\output_1D_athero_fluid.mat']);

param = declare_fixed_parameters_0D(ATHERO);

for i = 1:size(p_apv,2)
    Out = lumped_model(varpar(i),param);
    p_0D(i) = Out(3);
    p_1D(i) = p_apv(i);
end

ax2 = nexttile;hold all;
set(gca,'FontSize',14,'TickLabelInterpreter','latex');
title('atherosclerotic','Interpreter','latex');
text(-0.2,1,labels{2},'Units','normalized','FontSize',13,'fontweight','bold','FontName','times');
xlabel('$p_{\rm apv}$ (radial)','Interpreter','latex')
ylabel('$p_{\rm apv}$ (lumped-parameter)','Interpreter','latex')
[p_1Dsort,sortIdx] = sort(p_1D);
pdiffa = p_1Dsort-p_0D(sortIdx);
f5 = fill([-20 -15 15 15 10 -20],[-20 -20 10 15 15 -15],[0 0.8 0],'FaceAlpha',0.2,'EdgeColor','none');
f2 = fill([-20 -18 15 15 13 -20],[-20 -20 13 15 15 -18],[0 0.8 0],'FaceAlpha',0.45,'EdgeColor','none');
plot([-20 15],[-20 15],'linew',1.5,'color','k');
scatter(p_1Dsort,p_0D(sortIdx),10,'filled','+','MarkerEdgeColor',[0 0 0.8]);
text(-15,12,'{\it paired comparisons}','FontSize',13,'Interpreter','latex');
xlim([-20 15])
ylim([-20 15])
legend([f2 f5],{'$\pm$ 2 mmHg';'$\pm$ 5 mmHg'},'fontsize',14,'Interpreter','latex','location','southeast');


if EXPORT
    exportgraphics(gcf,[hompeath 'figures\papv_comparisons.png'],'Resolution',300);
end
