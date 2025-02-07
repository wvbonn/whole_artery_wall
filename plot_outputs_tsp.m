% plot max gradient magnitude and DC-CCL19 sensitivity distance

clear 
close all

currdir = pwd;
k = strfind(currdir,'\');
homepath = currdir(1:k(end));

ATHERO = 1==1;
EXPORT = 1==0;
labels = get_subplot_labels('a':'z',8);

param = declare_fixed_parameters(ATHERO);
if ATHERO
    % load numbers from the other config to get axis limits
    load('output_healthy_tsp.mat','Pe','Da_v','Da_l','varpar');
    Pe_2 = Pe; Da_v_2 = Da_v; Da_l_2 = Da_l; 
    for i = 1:numel(varpar); R_d_2(i) = varpar(i).R_d; end
    load('output_athero_tsp.mat','Pe','Da_v','Da_l','varpar','max_dc','max_dc_id','t_dc_ccl');
else
    % load numbers from the other config to get axis limits
    load('output_athero_tsp.mat','Pe','Da_v','Da_l','varpar');
    Pe_2 = Pe; Da_v_2 = Da_v; Da_l_2 = Da_l; 
    for i = 1:numel(varpar); R_d_2(i) = varpar(i).R_d; end
    load('output_healthy_tsp.mat','Pe','Da_v','Da_l','varpar','max_dc','max_dc_id','t_dc_ccl');
end
N = size(max_dc,1);
for i = 1:N; R_d(i) = varpar(i).R_d; end
n_D = size(max_dc,2);

xlims.Pe = [0.99*min([Pe Pe_2],[],'all') 1.01*max([Pe Pe_2],[],'all')];
xlims.Da_v = [0.99*min([Da_v Da_v_2],[],'all') 1.01*max([Da_v Da_v_2],[],'all')];
xlims.Da_l = [0.99*min([Da_l Da_l_2],[],'all') 1.01*max([Da_l Da_l_2],[],'all')];
ylims.dc = [0 13];
y_ticks = 0:2:10;
ylims.t_dC = [-0.4 1.2];

return;

%% plot max. absolute gradients against non-dimensional numbers

col = [[1 0.6 0.6];[0 0.6 1];[1 0.6 0.3]];
mkr = {'x';'o';'+';'s';'.'};
ls = {'-';'-';'--';'-.'};
txt = {'$D = 10^{-10}\,{\rm m}^2/{\rm s}$';
       '$D = 10^{-11}\,{\rm m}^2/{\rm s}$';
       '$D = 10^{-12}\,{\rm m}^2/{\rm s}$';
       'max. $||\nabla C||$ not at EEL'};
sz = 10; % size of scatter points

figure('position',[50 50 1000 800],'color','w');
tiledlayout(15,2);

% legend
nexttile([1 2]);hold all;
set(gca,'Visible','off');
set(gca,'XLim',[0 1]);
set(gca,'YLim',[0 1]);
x_legend = 0:0.235:1;
for i = 1:n_D
    scatter(x_legend(i),0.5,60,col(i,:),mkr{i},'linew',1.5);
    text(x_legend(i)+0.02,0.5,txt{i},'FontSize',15,'Interpreter','latex');
end
scatter(0.7,0.5,60,'k','.');
text(0.7+0.02,0.5,txt{4},'FontSize',16,'Interpreter','latex');

% Péclet
nexttile([7 1]);hold all;
set(gca,'fontsize',14,'TickLabelInterpreter','latex');
text(-0.15,1,labels{1},'Units','normalized','FontSize',13,'fontweight','bold','FontName','times');
xlabel('Pe','Interpreter','latex');
ylabel('${\rm max} ||\nabla C|| ({\rm mm}^{-1})$','Interpreter','latex');
ylim(ylims.dc)
set(gca,'XScale','log');
xl = xlims.Pe; xlim(xl);
xl_int = [ceil(log10(xl(1))) floor(log10(xl(end)))];
x_ticks = [];
for i = xl_int(1):xl_int(2)
    x_ticks = [x_ticks 10^i];
end
xticks(x_ticks);
yticks(y_ticks);
fill([1 1 xl(1) xl(1)],[ylims.dc fliplr(ylims.dc)],0.8*[1 1 1],'edgecolor','none','facealpha',0.2);
text(xl(1)+0.01*log10(diff(xl)),12,'diffusion','fontsize',13,'interpreter','latex');
text(xl(1)+0.5*log10(diff(xl)),12,'convection by interstitial fluid','fontsize',13,'interpreter','latex');
line([1 1],ylims.dc,'color','k','linew',1);
for i = 1:n_D
    scatter(Pe(:,i),max_dc(:,i),sz,col(i,:),mkr{i});
    scatter(Pe(max_dc_id(:,i)~=1,i),max_dc(max_dc_id(:,i)~=1,i),1/3*sz,'k','filled');
end

% vasa-vasorum Damköhler
ax_v = nexttile([7 1]);hold all;
set(gca,'fontsize',14,'TickLabelInterpreter','latex');
text(-0.15,1,labels{2},'Units','normalized','FontSize',13,'fontweight','bold','FontName','times');
xlabel('${\rm Da}_v$','Interpreter','latex');
ylabel('${\rm max} ||\nabla C|| ({\rm mm}^{-1})$','Interpreter','latex');
ylim(ylims.dc)
set(gca,'XScale','log');
xl = xlims.Da_v; xlim(xl);
xl_int = [ceil(log10(xl(1))) floor(log10(xl(end)))];
x_ticks = [];
for i = xl_int(1):xl_int(2)
    x_ticks = [x_ticks 10^i];
end
xticks(x_ticks);
yticks(y_ticks);
fill([1 1 xl(1) xl(1)],[ylims.dc fliplr(ylims.dc)],0.8*[1 1 1],'edgecolor','none','facealpha',0.2);
text(xl(1)+0.01*log10(diff(xl)),12,'diffusion','fontsize',13,'interpreter','latex');
text(xl(1)+0.5*log10(diff(xl)),12,'dilution by fluid from v.v.','fontsize',13,'interpreter','latex');
line([1 1],ylims.dc,'color','k','linew',1);
for i = 1:n_D
    scatter(Da_v(:,i),max_dc(:,i),sz,col(i,:),mkr{i});
    scatter(Da_v(max_dc_id(:,i)~=1,i),max_dc(max_dc_id(:,i)~=1,i),1/3*sz,'k','filled');
end

% lymphatic Damköhler
ax_l = nexttile([7 1]);hold all;
set(gca,'fontsize',14,'TickLabelInterpreter','latex');
text(-0.15,1,labels{3},'Units','normalized','FontSize',13,'fontweight','bold','FontName','times');
xlabel('${\rm Da}_{\ell}$','Interpreter','latex');
ylabel('${\rm max} ||\nabla C|| ({\rm mm}^{-1})$','Interpreter','latex');
ylim(ylims.dc)
set(gca,'XScale','log');
xl = xlims.Da_l; xlim(xl);
xl_int = [ceil(log10(xl(1))) floor(log10(xl(end)))];
x_ticks = [];
for i = xl_int(1):xl_int(2)
    x_ticks = [x_ticks 10^i];
end
xticks(x_ticks);
yticks(y_ticks);
fill([1 1 xl(1) xl(1)],[ylims.dc fliplr(ylims.dc)],0.8*[1 1 1],'edgecolor','none','facealpha',0.2);
text(xl(1)+0.01*log10(diff(xl)),12,'diffusion','fontsize',13,'interpreter','latex');
text(xl(1)+log10(diff(xl)),12,'convection into l.v.','fontsize',13,'interpreter','latex');
line([1 1],ylims.dc,'color','k','linew',1);
for i = 1:n_D
    scatter(Da_l(:,i),max_dc(:,i),sz,col(i,:),mkr{i});
    scatter(Da_l(max_dc_id(:,i)~=1,i),max_dc(max_dc_id(:,i)~=1,i),1/3*sz,'k','filled');
end


nexttile([7 1]);hold all;
set(gca,'fontsize',14,'TickLabelInterpreter','latex');
text(-0.15,1,labels{4},'Units','normalized','FontSize',13,'fontweight','bold','FontName','times');
xlabel('${\rm R}_{\rm d}$','Interpreter','latex');
ylabel('${\rm max} ||\nabla C|| ({\rm mm}^{-1})$','Interpreter','latex');
pos = get(gca,'position');
set(gca,'visible','off');
ax(1) = axes('position',[pos(1) pos(2) pos(3)/3 pos(4)]);hold(ax(1),'on');
ax(2) = axes('position',[pos(1)+pos(3)/3 pos(2) pos(3)/3 pos(4)]);hold(ax(2),'on');
ax(3) = axes('position',[pos(1)+2*pos(3)/3 pos(2) pos(3)/3 pos(4)]);hold(ax(3),'on');
set(ax(1),'XScale','log','YLim',ylims.dc); set(ax(1),'fontsize',15,'fontname','times');
set(ax(2),'XScale','log','YLim',ylims.dc); set(ax(2),'fontsize',15,'fontname','times');
set(ax(3),'XScale','log','YLim',ylims.dc); set(ax(3),'fontsize',15,'fontname','times');
ax(2).XLabel.String = '${\rm R}_{\rm d}$';
ax(2).XLabel.Interpreter = 'latex';
ax(1).YLabel.String = '${\rm max} ||\nabla C|| ({\rm mm}^{-1})$';
ax(1).YLabel.Interpreter = 'latex';
ax(1).YTick = y_ticks;
ax(2).YTick = [];
ax(3).YTick = [];
ax(1).XTickLabel = {'10^0';'10^1';'10^2|10^0'};
ax(2).XTickLabel = {'';'10^1';'10^2|10^0'};
ax(3).XTickLabel = {'';'';'10^2'};
for i = 1:n_D
    scatter(ax(i),R_d,max_dc(:,i),sz,col(i,:),mkr{i});
    scatter(ax(i),R_d(max_dc_id(:,i)~=1),max_dc(max_dc_id(:,i)~=1,i),1/3*sz,'k','filled');
end

ax_v.XLim = xlims.Da_v;
ax_l.XLim = xlims.Da_l;

if EXPORT
    if ATHERO
        exportgraphics(gcf,[homepath 'figures\max_dc_athero.png'],'Resolution',300);
    else
        exportgraphics(gcf,[homepath 'figures\max_dc_healthy.png'],'Resolution',300);
    end
end



%% plot dc-ccl19 sensitivity threshold distances

figure('position',[50 50 1000 800],'color','w');
tiledlayout(15,2);

nexttile([1 2]);hold all;
set(gca,'Visible','off');
set(gca,'XLim',[0 1]);
set(gca,'YLim',[0 1]);
x_legend = 0:1/3:1;
for i = 1:n_D
    scatter(x_legend(i),0.5,60,col(i,:),mkr{i},'linew',1.5);
    text(x_legend(i)+0.02,0.5,txt{i},'FontSize',15,'Interpreter','latex');
end

nexttile([7 1]);hold all;
set(gca,'fontsize',14,'TickLabelInterpreter','latex');
text(-0.15,1,labels{1},'Units','normalized','FontSize',13,'fontweight','bold','FontName','times');
xlabel('Pe','Interpreter','latex');
ylabel('$t^*_{\rm DC,PVAT}$','Interpreter','latex');
pos = get(gca,'position');
ylim(ylims.t_dC);
set(gca,'XScale','log');
xl = xlims.Pe; xlim(xl);
line(xl,[0 0],'color','k','lines','--','linew',1.5);
annotation('textbox',[pos(1) pos(2)+0.14*pos(4) pos(3) pos(4)],'String','adventitia','FontSize',14,'FontName','times',...
    'EdgeColor','none','VerticalAlignment','bottom','horizontalalignment','left','FitBoxToText','on');
annotation('textbox',[pos(1) pos(2)+0.25*pos(4) pos(3) pos(4)],'String','PVAT','FontSize',14,'FontName','times',...
    'EdgeColor','none','VerticalAlignment','bottom','horizontalalignment','left','FitBoxToText','on');
xl_int = [ceil(log10(xl(1))) floor(log10(xl(end)))];
x_ticks = [];
for i = xl_int(1):xl_int(2)
    x_ticks = [x_ticks 10^i];
end
xticks(x_ticks);
fill([1 1 xl(1) xl(1)],[ylims.t_dC fliplr(ylims.t_dC)],0.8*[1 1 1],'edgecolor','none','facealpha',0.2);
text(xl(1)+0.01*log10(diff(xl)),1.1,'diffusion','fontsize',13,'interpreter','latex');
text(xl(1)+0.5*log10(diff(xl)),1.1,'convection by interstitial fluid','fontsize',13,'interpreter','latex');
line([1 1],ylims.t_dC,'color','k','linew',1);
for i = 1:n_D
    scatter(Pe(:,i),t_dc_ccl(:,i),sz,col(i,:),mkr{i});
end


ax_v = nexttile([7 1]);hold all;
set(gca,'fontsize',14,'TickLabelInterpreter','latex');
text(-0.15,1,labels{2},'Units','normalized','FontSize',13,'fontweight','bold','FontName','times');
xlabel('${\rm Da}_v$','Interpreter','latex');
ylabel('$t^*_{\rm DC,PVAT}$','Interpreter','latex');
pos = get(ax_v,'position');
ylim(ylims.t_dC);
set(gca,'XScale','log');
xl = xlims.Da_v; xlim(xl);
line(xl,[0 0],'color','k','lines','--','linew',1.5);
annotation('textbox',[pos(1) pos(2)+0.14*pos(4) pos(3) pos(4)],'String','adventitia','FontSize',14,'FontName','times',...
    'EdgeColor','none','VerticalAlignment','bottom','horizontalalignment','left','FitBoxToText','on');
annotation('textbox',[pos(1) pos(2)+0.25*pos(4) pos(3) pos(4)],'String','PVAT','FontSize',14,'FontName','times',...
    'EdgeColor','none','VerticalAlignment','bottom','horizontalalignment','left','FitBoxToText','on');
xl_int = [ceil(log10(xl(1))) floor(log10(xl(end)))];
x_ticks = [];
for i = xl_int(1):xl_int(2)
    x_ticks = [x_ticks 10^i];
end
xticks(x_ticks);
fill([1 1 xl(1) xl(1)],[ylims.t_dC fliplr(ylims.t_dC)],0.8*[1 1 1],'edgecolor','none','facealpha',0.2);
text(xl(1)+0.01*log10(diff(xl)),1.1,'diffusion','fontsize',13,'interpreter','latex');
text(xl(1)+0.5*log10(diff(xl)),1.1,'dilution by fluid from v.v.','fontsize',13,'interpreter','latex');
line([1 1],ylims.t_dC,'color','k','linew',1);
for i = 1:n_D
    scatter(Da_v(:,i),t_dc_ccl(:,i),sz,col(i,:),mkr{i});
end


ax_l = nexttile([7 1]);hold all;
set(gca,'fontsize',14,'TickLabelInterpreter','latex');
text(-0.15,1,labels{3},'Units','normalized','FontSize',13,'fontweight','bold','FontName','times');
xlabel('${\rm Da}_{\ell}$','Interpreter','latex');
ylabel('$t^*_{\rm DC,PVAT}$','Interpreter','latex');
pos = get(ax_l,'position');
ylim(ylims.t_dC);
set(gca,'XScale','log');
xl = xlims.Da_l; xlim(xl);
line(xl,[0 0],'color','k','lines','--','linew',1.5);
annotation('textbox',[pos(1) pos(2)+0.14*pos(4) pos(3) pos(4)],'String','adventitia','FontSize',14,'FontName','times',...
    'EdgeColor','none','VerticalAlignment','bottom','horizontalalignment','left','FitBoxToText','on');
annotation('textbox',[pos(1) pos(2)+0.25*pos(4) pos(3) pos(4)],'String','PVAT','FontSize',14,'FontName','times',...
    'EdgeColor','none','VerticalAlignment','bottom','horizontalalignment','left','FitBoxToText','on');
xl_int = [ceil(log10(xl(1))) floor(log10(xl(end)))];
x_ticks = [];
for i = xl_int(1):xl_int(2)
    x_ticks = [x_ticks 10^i];
end
xticks(x_ticks);
fill([1 1 xl(1) xl(1)],[ylims.t_dC fliplr(ylims.t_dC)],0.8*[1 1 1],'edgecolor','none','facealpha',0.2);
text(xl(1)+0.01*log10(diff(xl)),1.1,'diffusion','fontsize',13,'interpreter','latex');
text(xl(1)+log10(diff(xl)),1.1,'convection into l.v.','fontsize',13,'interpreter','latex');
line([1 1],ylims.t_dC,'color','k','linew',1);
for i = 1:n_D
    scatter(Da_l(:,i),t_dc_ccl(:,i),sz,col(i,:),mkr{i});
end

ax0 = nexttile([7 1]);hold all;
set(ax0,'fontsize',15,'fontname','times');
text(-0.15,1,labels{4},'Units','normalized','FontSize',13,'fontweight','bold','FontName','times');
set(ax0,'XTick',[],'YTick',[])
xlabel('${\rm R}_{\rm d}$','Interpreter','latex');
ylabel('$t^*_{\rm DC,PVAT}$','Interpreter','latex');
pos = get(ax0,'position');
set(ax0,'visible','off');
ax(1) = axes('position',[pos(1) pos(2) pos(3)/3 pos(4)]);hold(ax(1),'on');
ax(2) = axes('position',[pos(1)+pos(3)/3 pos(2) pos(3)/3 pos(4)]);hold(ax(2),'on');
ax(3) = axes('position',[pos(1)+2*pos(3)/3 pos(2) pos(3)/3 pos(4)]);hold(ax(3),'on');
set(ax(1),'XScale','log','YLim',ylims.t_dC); set(ax(1),'fontsize',15,'fontname','times');
set(ax(2),'XScale','log','YLim',ylims.t_dC); set(ax(2),'fontsize',15,'fontname','times');
set(ax(3),'XScale','log','YLim',ylims.t_dC); set(ax(3),'fontsize',15,'fontname','times');
ax(2).XLabel.String = '${\rm R}_{\rm d}$';
ax(2).XLabel.Interpreter = 'latex';
ax(1).YLabel.String = '$t^*_{\rm DC,PVAT}$';
ax(1).YLabel.Interpreter = 'latex';
ax(2).YTick = [];
ax(3).YTick = [];
ax(1).XTickLabel = {'10^0';'10^1';'10^2|10^0'};
ax(2).XTickLabel = {'';'10^1';'10^2|10^0'};
ax(3).XTickLabel = {'';'';'10^2'};
for i = 1:n_D
    scatter(ax(i),R_d,t_dc_ccl(:,i),sz,col(i,:),mkr{i});
end
line(ax(1),ax(1).XLim,[0 0],'color','k','lines','--','linew',1.5);
line(ax(2),ax(2).XLim,[0 0],'color','k','lines','--','linew',1.5);
line(ax(3),ax(3).XLim,[0 0],'color','k','lines','--','linew',1.5);
annotation('textbox',[pos(1) pos(2)+0.14*pos(4) pos(3) pos(4)],'String','adventitia','FontSize',14,'FontName','times',...
    'EdgeColor','none','VerticalAlignment','bottom','horizontalalignment','left','FitBoxToText','on');
annotation('textbox',[pos(1) pos(2)+0.25*pos(4) pos(3) pos(4)],'String','PVAT','FontSize',14,'FontName','times',...
    'EdgeColor','none','VerticalAlignment','bottom','horizontalalignment','left','FitBoxToText','on');

ax_v.XLim = xlims.Da_v;
ax_l.XLim = xlims.Da_l;

if EXPORT
    if ATHERO
        exportgraphics(gcf,[homepath 'figures\t_dc_pvat_athero.png'],'Resolution',300);
    else
        exportgraphics(gcf,[homepath 'figures\t_dc_pvat_healthy.png'],'Resolution',300);
    end
end

