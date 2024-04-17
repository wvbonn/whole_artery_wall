clear
close all

homepath = 'C:\Users\Willy\Work\PhD\ATLO\clean\';
EXPORT = 1==1;

AR = 1.5; % aspect ratio of figure (needs to be AR*screen size in height)
% axis dimensions
axd.x0 = [0.1 0.54];
axd.w = 0.36;
axd.y0 = [0.11 0.38 0.65 0.89];
axd.h = 0.16;
% declare plot styles
l.len = 0.025; % length of legend line
% ADAPT COLOR TO FIGS 13-16: [1 0.6 0.6];[0 0.6 1];[0.6 0.6 0]
col = [0.3;0.7].*ones(2,3); % color, k_im7
ls = {':';'-'}; % line style, R_d
ms = {'*';'^'}; % marker style, q_lv

% permutation of parameters (verify varpar in the data output file if needed)
col_id = [1 2 1 2 1 2 1 2]; % k_im
ls_id = [1 1 2 2 1 1 2 2]; % R_d
ms_id = [1 1 1 1 2 2 2 2]; % q_lv

ATHERO = [1==0,1==1];
% load and plot data
r_eel = 2;

r_lu(1) = r_eel-0.34;
r_apv(1) = r_eel+0.54;
r_ext(1) = r_apv(1)+0.88;
k_im(1,:) = [0.6 1.5];
q_lv(1,:) = [0.6 1.4];

r_lu(2) = r_eel-0.78;
r_apv(2) = r_eel+0.92;
r_ext(2) = r_apv(2)+0.88;
k_im(2,:) = [1.2 2.8];
q_lv(2,:) = [0.2 1.2];
    
%% prepare figure
figure('position',[50 0 1200/AR 1000],'color','w');
hold on;
set(gca,'visible','off');
% tiledlayout(21,2,"TileSpacing","compact");
labels = get_subplot_labels('a':'z',6);

axl = axes('Position',[axd.x0(1) axd.y0(4) 1-axd.x0(1) 2*axd.h/3]);hold on;
xlim([0 1])
ylim([0 1])
set(gca,'visible','off');
text(0,0.9,'Color:','FontSize',1/AR*15,'Interpreter','latex');
line([0 l.len],0.67*[1 1],'color',col(2,:),'linew',1.5);
text(0.26,0.75,{['$k_{\rm im} = ' num2str(k_im(1,2)) '\,10^{-18}\,{\rm m}^2$ (healthy)'];['$' num2str(k_im(2,2)) '\,10^{-18}\,{\rm m}^2$ (athero.)']},...
    'FontSize',1/AR*15,'Interpreter','latex','HorizontalAlignment','right','VerticalAlignment','top');
line([0 l.len],0.17*[1 1],'color',col(1,:),'linew',1.5);
text(0.26,0.25,{['$k_{\rm im} = ' num2str(k_im(1,1)) '\,10^{-18}\,{\rm m}^2$ (healthy)'];['$' num2str(k_im(2,1)) '\,10^{-18}\,{\rm m}^2$ (athero.)']},...
    'FontSize',1/AR*15,'Interpreter','latex','HorizontalAlignment','right','VerticalAlignment','top');
text(0.34,0.9,'Line style:','FontSize',1/AR*15,'Interpreter','latex');
line(0.34+[0 l.len],0.67*[1 1],'color',mean(col,1),'linew',1.5,'lines',ls{2});
text(0.38,0.75,'${\rm R}_{\rm d}$ = 30','FontSize',1/AR*15,'Interpreter','latex','VerticalAlignment','top');
line(0.34+[0 l.len],0.17*[1 1],'color',mean(col,1),'linew',1.5,'lines',ls{1});
text(0.38,0.25,'${\rm R}_{\rm d}$ = 3','FontSize',1/AR*15,'Interpreter','latex','VerticalAlignment','top');
text(0.54,0.9,'Marker:','FontSize',1/AR*15,'Interpreter','latex');
plot(0.54+0.025,0.67,'Marker',ms{2},'color',mean(col,1),'MarkerSize',1/AR*10,'linew',1);
text(0.82,0.75,{['$\tilde{q}_{\rm lv} = ' num2str(q_lv(1,2)) '\,10^{-11}\,{\rm m}^2$/s (healthy)'];['$' num2str(q_lv(2,2)) '\,10^{-11}\,{\rm m}^2$/s (athero.)']},...
    'FontSize',1/AR*15,'Interpreter','latex','HorizontalAlignment','right','VerticalAlignment','top');
plot(0.54+0.025,0.17,'Marker',ms{1},'color',mean(col,1),'MarkerSize',1/AR*10,'linew',1);
text(0.82,0.25,{['$\tilde{q}_{\rm lv} = ' num2str(q_lv(1,1)) '\,10^{-11}\,{\rm m}^2$/s (healthy)'];['$' num2str(q_lv(2,1)) '\,10^{-11}\,{\rm m}^2$/s (athero.)']},...
    'FontSize',1/AR*15,'Interpreter','latex','HorizontalAlignment','right','VerticalAlignment','top');


ax(1,1) = axes('Position',[axd.x0(1) axd.y0(3) axd.w axd.h]);hold on;
set(gca,'fontsize',1/AR*15,'TickLabelInterpreter','latex');
xlabel('r (mm)','Interpreter','latex');
ylabel('c (-)','Interpreter','latex');
xlim([r_lu(1) r_ext(1)])
ylim([0 1])
title('$D = 10^{-10}\,{\rm m}^2$/s (healthy)','FontSize',1/AR*15,'Interpreter','latex');
text(-0.12,1,labels{1},'Units','normalized','FontSize',1/AR*15,'fontweight','bold','fontname','times');

ax(1,2) = axes('Position',[axd.x0(2) axd.y0(3) axd.w axd.h]);hold on;
set(gca,'fontsize',1/AR*15,'TickLabelInterpreter','latex');
xlabel('r (mm)','Interpreter','latex');
ylabel('c (-)','Interpreter','latex');
xlim([r_lu(2) r_ext(2)])
ylim([0 1])
title('$D = 10^{-10}\,{\rm m}^2$/s (athero.)','FontSize',1/AR*15,'Interpreter','latex');
text(-0.12,1,labels{2},'Units','normalized','FontSize',1/AR*15,'fontweight','bold','fontname','times');

ax(2,1) = axes('Position',[axd.x0(1) axd.y0(2) axd.w axd.h]);hold on;
set(gca,'fontsize',1/AR*15,'TickLabelInterpreter','latex');
xlabel('r (mm)','Interpreter','latex');
ylabel('c (-)','Interpreter','latex');
xlim([r_lu(1) r_ext(1)])
ylim([0 1])
title('$D = 10^{-11}\,{\rm m}^2$/s (healthy)','FontSize',1/AR*15,'Interpreter','latex');
text(-0.12,1,labels{3},'Units','normalized','FontSize',1/AR*15,'fontweight','bold','fontname','times');

ax(2,2) = axes('Position',[axd.x0(2) axd.y0(2) axd.w axd.h]);hold on;
set(gca,'fontsize',1/AR*15,'TickLabelInterpreter','latex');
xlabel('r (mm)','Interpreter','latex');
ylabel('c (-)','Interpreter','latex');
xlim([r_lu(2) r_ext(2)])
ylim([0 1])
title('$D = 10^{-11}\,{\rm m}^2$/s (athero.)','FontSize',1/AR*15,'Interpreter','latex');
%text(2.8,0.9,{'D = 10^{-11} m^2/s'},'FontSize',1/AR*15,'fontname','times')
text(-0.12,1,labels{4},'Units','normalized','FontSize',1/AR*15,'fontweight','bold','fontname','times');

ax(3,1) = axes('Position',[axd.x0(1) axd.y0(1) axd.w axd.h]);hold on;
set(gca,'fontsize',1/AR*15,'TickLabelInterpreter','latex');
xlabel('r (mm)','Interpreter','latex');
ylabel('c (-)','Interpreter','latex');
xlim([r_lu(1) r_ext(1)])
ylim([0 1.1])
title('$D = 10^{-12}\,{\rm m}^2$/s (healthy)','FontSize',1/AR*15,'Interpreter','latex');
text(-0.12,1,labels{5},'Units','normalized','FontSize',1/AR*15,'fontweight','bold','fontname','times');

ax(3,2) = axes('Position',[axd.x0(2) axd.y0(1) axd.w axd.h]);hold on;
set(gca,'fontsize',1/AR*15,'TickLabelInterpreter','latex');
xlabel('r (mm)','Interpreter','latex');
ylabel('c (-)','Interpreter','latex');
xlim([r_lu(2) r_ext(2)])
ylim([0 1.1])
title('$D = 10^{-12}\,{\rm m}^2$/s (athero.)','FontSize',1/AR*15,'Interpreter','latex');
text(-0.12,1,labels{6},'Units','normalized','FontSize',1/AR*15,'fontweight','bold','fontname','times');

%% load data

load('nonrnd_1D_healthy_tsp.mat');
for i = 1:numel(r)
    plot(ax(1,1),r{i},c{i,1}/c{i,1}(1),'color',col(col_id(i),:),'lines',ls{ls_id(i)},'linew',1);
    scatter(ax(1,1),r{i}(3*i:10:end),c{i,1}(3*i:10:end)/c{i,1}(1),1/AR*30,col(col_id(i),:),'marker',ms{ms_id(i)},'linew',1);
    plot(ax(2,1),r{i},c{i,2}/c{i,2}(1),'color',col(col_id(i),:),'lines',ls{ls_id(i)},'linew',1);
    scatter(ax(2,1),r{i}(3*i:10:end),c{i,2}(3*i:10:end)/c{i,2}(1),1/AR*30,col(col_id(i),:),'marker',ms{ms_id(i)},'linew',1);
    plot(ax(3,1),r{i},c{i,3}/c{i,3}(1),'color',col(col_id(i),:),'lines',ls{ls_id(i)},'linew',1);
    scatter(ax(3,1),r{i}(3*i:10:end),c{i,3}(3*i:10:end)/c{i,3}(1),1/AR*30,col(col_id(i),:),'marker',ms{ms_id(i)},'linew',1);
end

load('nonrnd_1D_athero_tsp.mat');
for i = 1:numel(r)
    plot(ax(1,2),r{i},c{i,1}/c{i,1}(1),'color',col(col_id(i),:),'lines',ls{ls_id(i)},'linew',1);
    scatter(ax(1,2),r{i}(3*i:10:end),c{i,1}(3*i:10:end)/c{i,1}(1),1/AR*30,col(col_id(i),:),'marker',ms{ms_id(i)},'linew',1);
    plot(ax(2,2),r{i},c{i,2}/c{i,2}(1),'color',col(col_id(i),:),'lines',ls{ls_id(i)},'linew',1);
    scatter(ax(2,2),r{i}(3*i:10:end),c{i,2}(3*i:10:end)/c{i,2}(1),1/AR*30,col(col_id(i),:),'marker',ms{ms_id(i)},'linew',1);
    plot(ax(3,2),r{i},c{i,3}/c{i,3}(1),'color',col(col_id(i),:),'lines',ls{ls_id(i)},'linew',1);
    scatter(ax(3,2),r{i}(3*i:10:end),c{i,3}(3*i:10:end)/c{i,3}(1),1/AR*30,col(col_id(i),:),'marker',ms{ms_id(i)},'linew',1);
end

% draw domain boundaries
str_a = {'adventitia';'adventitia'};
for j = 1:2
    xl = get(ax(1,j),'XLim');
    for i = 1:3
        yl = get(ax(i,j),'YLim');
        if i==1; text(ax(i,j),(xl(1)+r_eel)/2,yl(1)+0.95*diff(yl),{'inner';'layers'},'fontsize',1/AR*15,'Interpreter','latex','HorizontalAlignment','center','VerticalAlignment','top'); end
        line(ax(i,j),r_eel*[1 1],yl,'color','k');
        if i==1; text(ax(i,j),(r_eel+r_apv(j))/2,yl(1)+1*diff(yl),str_a{j},'fontsize',1/AR*15,'Interpreter','latex','HorizontalAlignment','center','VerticalAlignment','top'); end
        line(ax(i,j),r_apv(j)*[1 1],yl,'color','k');
        %text(ax1,0.8*r_apv+0.2*xl(2),yl(1)+0.04*diff(yl),'PVAT','fontsize',1/AR*15,'Interpreter','latex','HorizontalAlignment','center','VerticalAlignment','bottom');
    end
end

if EXPORT
    exportgraphics(gcf,[homepath 'figures\nonrnd_tsp_variables.png'],'Resolution',300);
end