clear
close all

homepath = 'C:\Users\Willy\Work\PhD\ATLO\clean\';
EXPORT = 1==1;
% declare plot styles
col = [0.7;0.3;0.5].*ones(3,3); % color, k_im
ls = {'-';':'}; % line style, k_a
ms = {'*';'^'}; % marker style, q_lv

% permutation of parameters (verify varpar in the data output file if needed)
col_id = [1 2 1 2 1 2 1 2]; % k_im
ls_id = [1 1 2 2 1 1 2 2]; % k_a
ms_id = [1 1 1 1 2 2 2 2]; % q_lv

% load and plot data
load('nonrnd_1D_healthy_fluid.mat');

r_eel = 2;
r_lu = r_eel-0.34;
r_apv = r_eel+0.56;
r_ext = r_apv+0.9;

figure('position',[50 0 1200 800],'color','w');
tiledlayout(2,2);
labels = get_subplot_labels('a':'z',6);

ax1 = nexttile;hold all;
set(gca,'fontsize',15,'TickLabelInterpreter','latex');
xlabel('$r$ (mm)','Interpreter','latex');
ylabel('$p$ (mmHg)','Interpreter','latex');
xlim([r_lu r_ext])
ylim([-20 20]);
text(-0.15,1,labels{1},'Units','normalized','FontSize',14,'fontweight','bold','fontname','times');

ax2 = nexttile;hold all;
set(gca,'fontsize',15,'TickLabelInterpreter','latex');
xlabel('$r$ (mm)','Interpreter','latex');
ylabel('$v$ ($\mu$m/s)','Interpreter','latex');
xlim([r_lu r_ext])
ylim([0 0.045])
text(-0.15,1,labels{2},'Units','normalized','FontSize',14,'fontweight','bold','fontname','times');

ax3 = nexttile;hold all;
set(gca,'fontsize',15,'TickLabelInterpreter','latex');
xlabel('$r$ (mm)','Interpreter','latex');
ylabel('$q_{\rm vv}$ (${\rm s}^{-1}$)','Interpreter','latex');
xlim([r_lu r_ext])
text(-0.15,1,labels{3},'Units','normalized','FontSize',14,'fontweight','bold','fontname','times');

axl = nexttile;hold all;
set(gca,'Visible','off');
xlim([0 1])
ylim([0 1])
text(0,0.9,'Color:','FontSize',15,'Interpreter','latex');
line([0 0.1],0.8*[1 1],'color',col(1,:),'linew',1.5);
text(0.15,0.8,'$k_{\rm im} = 0.7\,10^{-18}\,{\rm m}^2$','FontSize',15,'Interpreter','latex');
line([0 0.1],0.7*[1 1],'color',col(2,:),'linew',1.5);
text(0.15,0.7,'$k_{\rm im} = 1.2\,10^{-18}\,{\rm m}^2$','FontSize',15,'Interpreter','latex');
text(0,0.6,'Line style:','FontSize',15,'Interpreter','latex');
line([0 0.1],0.5*[1 1],'color',col(3,:),'linew',1.5,'lines',ls{1});
text(0.15,0.5,'$k_{\rm a} = 5\,10^{-18}\,{\rm m}^2$','FontSize',15,'Interpreter','latex');
line([0 0.1],0.4*[1 1],'color',col(3,:),'linew',1.5,'lines',ls{2});
text(0.15,0.4,'$k_{\rm a} = 10^{-16}\,{\rm m}^2$','FontSize',15,'Interpreter','latex');
text(0,0.3,'Marker:','FontSize',15,'Interpreter','latex');
plot(0.05,0.2,'Marker',ms{1},'color',col(3,:),'MarkerSize',8,'linew',1.5);
text(0.15,0.2,'$\tilde{q}_{\rm lv} = 5\,10^{-12}\,{\rm m}^2/{\rm s}$','FontSize',15,'Interpreter','latex');
plot(0.05,0.1,'Marker',ms{2},'color',col(3,:),'MarkerSize',8,'linew',1.5);
text(0.15,0.1,'$\tilde{q}_{\rm lv} = 7\,10^{-12}\,{\rm m}^2/{\rm s}$','FontSize',15,'Interpreter','latex');


for i = 1:numel(r)
    plot(ax1,r{i},p{i},'color',col(col_id(i),:),'lines',ls{ls_id(i)},'linew',1.5);
    scatter(ax1,r{i}(2*i:10:end),p{i}(2*i:10:end),30,col(col_id(i),:),'marker',ms{ms_id(i)},'linew',1.5);
    plot(ax2,r{i},v{i},'color',col(col_id(i),:),'lines',ls{ls_id(i)},'linew',1.5);
    scatter(ax2,r{i}(2*i:10:end),v{i}(2*i:10:end),30,col(col_id(i),:),'marker',ms{ms_id(i)},'linew',1.5);
    id_eel = find(r{i}>2,1);
    plot(ax3,r{i}(id_eel:end),Q_vv{i}(id_eel:end),'color',col(col_id(i),:),'lines',ls{ls_id(i)},'linew',1.5);
    scatter(ax3,r{i}(id_eel+2*i:10:end),Q_vv{i}(id_eel+2*i:10:end),30,col(col_id(i),:),'marker',ms{ms_id(i)},'linew',1.5);
end

for j = [1 5]
    plot(ax3,r{j}(id_eel:end),Q_lv{j}(id_eel:end),'color',[0 0.7 0],'lines','-','linew',1);
    scatter(ax3,r{j}(id_eel+j:10:end),Q_lv{j}(id_eel+2*j:10:end),30,[0 0.7 0],'marker',ms{ms_id(j)},'linew',1.5);
end
yl = get(ax3,'YLim');
text(ax3,r{1}(id_eel),Q_lv{1}(id_eel)+0.07*diff(yl),'$-q_{\rm lv}$','FontSize',15,'Interpreter','latex');
text(ax3,r{5}(id_eel),Q_lv{5}(id_eel)+0.07*diff(yl),'$-q_{\rm lv}$','FontSize',15,'Interpreter','latex');

% draw domain boundaries
xl = get(ax1,'XLim');

yl = get(ax1,'YLim');
text(ax1,(xl(1)+r_eel)/2,yl(1)+0.04*diff(yl),{'inner';'layers';'($p_{\rm lu}$ =';'90 mmHg)'},'fontsize',13,'interpreter','latex','HorizontalAlignment','center','VerticalAlignment','bottom');
line(ax1,r_eel*[1 1],yl,'color','k');
text(ax1,(r_eel+r_apv)/2,yl(1)+0.04*diff(yl),'adventitia','fontsize',13,'interpreter','latex','HorizontalAlignment','center','VerticalAlignment','bottom');
line(ax1,r_apv*[1 1],yl,'color','k');
text(ax1,0.8*r_apv+0.2*xl(2),yl(1)+0.04*diff(yl),'PVAT','fontsize',13,'interpreter','latex','HorizontalAlignment','center','VerticalAlignment','bottom');

yl = get(ax2,'YLim');
text(ax2,(xl(1)+r_eel)/2,yl(1)+0.04*diff(yl),{'inner';'layers'},'fontsize',13,'interpreter','latex','HorizontalAlignment','center','VerticalAlignment','bottom');
line(ax2,r_eel*[1 1],yl,'color','k');
text(ax2,(r_eel+r_apv)/2,yl(1)+0.04*diff(yl),'adventitia','fontsize',13,'interpreter','latex','HorizontalAlignment','center','VerticalAlignment','bottom');
line(ax2,r_apv*[1 1],yl,'color','k');
text(ax2,0.8*r_apv+0.2*xl(2),yl(1)+0.04*diff(yl),'PVAT','fontsize',13,'interpreter','latex','HorizontalAlignment','center','VerticalAlignment','bottom');

yl = get(ax3,'YLim');
text(ax3,(xl(1)+r_eel)/2,yl(1)+0.04*diff(yl),{'inner';'layers'},'fontsize',13,'interpreter','latex','HorizontalAlignment','center','VerticalAlignment','bottom');
line(ax3,r_eel*[1 1],yl,'color','k');
text(ax3,(r_eel+r_apv)/2,yl(1)+0.04*diff(yl),'adventitia','fontsize',13,'interpreter','latex','HorizontalAlignment','center','VerticalAlignment','bottom');
line(ax3,r_apv*[1 1],yl,'color','k');
text(ax3,0.8*r_apv+0.2*xl(2),yl(1)+0.04*diff(yl),'PVAT','fontsize',13,'interpreter','latex','HorizontalAlignment','center','VerticalAlignment','bottom');

if EXPORT
    exportgraphics(gcf,[homepath 'figures\nonrnd_fluid_variables.png'],'Resolution',300);
end