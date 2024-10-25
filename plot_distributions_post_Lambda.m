clear
close all

currdir = pwd;
k = strfind(currdir,'\');
homepath = currdir(1:k(end));

EXPORT = 1==0;
L_lim = [0.9 1.6];

load('unfiltered_distributions.mat');
k_i_0 = k_i; l_pv_0 = l_pv; n_v_0 = n_v; n_l_0 = n_l; q_l_0 = q_l;

load(['filtered_Lambda_' num2str(L_lim(1)) '_' num2str(L_lim(2)) '.mat']);

n = 10000; % number of draws
nBins = 50;
col.Hbase = [0 0.6 0];
col.HL = col.Hbase;
ls.L = '--';
col.Abase = [1 0 0.7];
col.AL = col.Abase;

figure('position',[50 50 800 650],'color','w');
tiledlayout(11,6);
labels = get_subplot_labels('a':'z',6);

%% legend

ax0 = nexttile([1 6]);hold all;
set(gca,'visible','off');
xlim([0 1]);
ylim([0 1]);
line([0 0.08],0*[1 1],'color',col.Hbase,'linew',2);
text(0.1,0,'healthy','VerticalAlignment','middle','FontSize',14,'Interpreter','latex');
plot([0.27 0.31 0.35],0*[1 1 1],'color',col.Abase,'linew',2,'marker','*','markerindices',2);
text(0.37,0,'atherosclerotic','VerticalAlignment','middle','FontSize',14,'Interpreter','latex');
line([0.68 0.8],0*[1 1],'color',0.4*[1 1 1],'linew',2,'lines',ls.L);
text(0.82,0,'physiologic $\Lambda$','VerticalAlignment','middle','FontSize',14,'Interpreter','latex');
% line([0.75 0.8],0.5*[1 1],'color',[0.8 0 0],'linew',1.5);
% text(0.82,0.5,'physiol. interval','VerticalAlignment','middle','FontSize',12);


%% v.v. conductivity distributions
ax3 = nexttile([5 3]);hold all;
set(gca,'fontsize',15,'TickLabelInterpreter','latex');
text(-0.1,1,labels{1},'Units','normalized','FontSize',14,'fontweight','bold','fontname','times');
xlabel('$l_{p,v}$ (cm/s/mmHg)','Interpreter','latex');
set(gca,'XLim',[-7.3 -6.4]);
xticks([log10(6e-8) -7 log10(2e-7)]);
xticklabels({'$6\,10^{-8}$';'$10^{-7}$';'$2\,10^{-7}$'});
ylabel('pdf','Interpreter','latex');
yticks([])

l_pv_H = random(l_pv_0.base,[1 n]);
l_pv_A = random(l_pv_0.atherofact,[1 n]).*l_pv_H;
h = histogram(log10(l_pv_H),nBins,'FaceColor',col.Hbase,'BinLimits',[min(log10(l_pv_H(:))) max(log10(l_pv_H(:)))],'Normalization','pdf','Visible','off');
plot((h.BinEdges(1:end-1)+h.BinEdges(2:end))/2,smooth(h.Values,5),'linew',2,'color',col.Hbase);
h = histogram(log10(l_pv_A),nBins,'FaceColor',col.Abase,'BinLimits',[min(log10(l_pv_A(:))) max(log10(l_pv_A(:)))],'Normalization','pdf','Visible','off');
plot((h.BinEdges(1:end-1)+h.BinEdges(2:end))/2,smooth(h.Values,5),'linew',2,'color',col.Abase,'marker','*','markerindices',5:5:nBins);

l_pv_L_H = random(l_pv.base,[1 n]);
l_pv_L_A = random(l_pv_0.atherofact,[1 n]).*l_pv_L_H;
h = histogram(log10(l_pv_L_H),nBins,'FaceColor',col.Hbase,'BinLimits',[min(log10(l_pv_L_H(:))) max(log10(l_pv_L_H(:)))],'Normalization','pdf','Visible','off');
plot((h.BinEdges(1:end-1)+h.BinEdges(2:end))/2,smooth(h.Values,5),'linew',2,'color',col.Hbase,'lines',ls.L);
h = histogram(log10(l_pv_L_A),nBins,'FaceColor',col.Abase,'BinLimits',[min(log10(l_pv_L_A(:))) max(log10(l_pv_L_A(:)))],'Normalization','pdf','Visible','off');
plot((h.BinEdges(1:end-1)+h.BinEdges(2:end))/2,smooth(h.Values,5),'linew',2,'color',col.Abase,'lines',ls.L,'marker','*','markerindices',2:5:nBins);


%% v.v. density
ax4 = nexttile([5 3]);hold all;
set(gca,'fontsize',15,'TickLabelInterpreter','latex');
text(-0.1,1,labels{2},'Units','normalized','FontSize',14,'fontweight','bold','fontname','times');
xlabel('$N_v ({\rm mm}^{-2})$','Interpreter','latex');
xticks([1 log10(20) log10(50) 2]);
xticklabels({'$10^1$';'$2\,10^1$';'$5\,10^1$';'$10^2$'});
ylabel('pdf','Interpreter','latex');
yticks([])

n_v_H = random(n_v_0.base,[1 n]);
n_v_A = random(n_v_0.atherofact,[1 n]).*n_v_H;
h = histogram(log10(n_v_H),nBins,'FaceColor',col.Hbase,'BinLimits',[min(log10(n_v_H(:))) max(log10(n_v_H(:)))],'Normalization','pdf','Visible','off');
plot((h.BinEdges(1:end-1)+h.BinEdges(2:end))/2,smooth(h.Values,5),'linew',2,'color',col.Hbase);
h = histogram(log10(n_v_A),nBins,'FaceColor',col.Abase,'BinLimits',[min(log10(n_v_A(:))) max(log10(n_v_A(:)))],'Normalization','pdf','Visible','off');
plot((h.BinEdges(1:end-1)+h.BinEdges(2:end))/2,smooth(h.Values,5),'linew',2,'color',col.Abase,'marker','*','markerindices',5:5:nBins);

n_v_L_H = random(n_v.base,[1 n]);
n_v_L_A = random(n_v_0.atherofact,[1 n]).*n_v_L_H;
h = histogram(log10(n_v_L_H),nBins,'FaceColor',col.Hbase,'BinLimits',[min(log10(n_v_L_H(:))) max(log10(n_v_L_H(:)))],'Normalization','pdf','Visible','off');
plot((h.BinEdges(1:end-1)+h.BinEdges(2:end))/2,smooth(h.Values,5),'linew',2,'color',col.Hbase,'lines',ls.L);
h = histogram(log10(n_v_L_A),nBins,'FaceColor',col.Abase,'BinLimits',[min(log10(n_v_L_A(:))) max(log10(n_v_L_A(:)))],'Normalization','pdf','Visible','off');
plot((h.BinEdges(1:end-1)+h.BinEdges(2:end))/2,smooth(h.Values,5),'linew',2,'color',col.Abase,'lines',ls.L,'marker','*','markerindices',5:5:nBins);


%% l.v. density
ax5 = nexttile([5 3]);hold all;
set(gca,'fontsize',15,'TickLabelInterpreter','latex');
text(-0.1,1,labels{3},'Units','normalized','FontSize',14,'fontweight','bold','fontname','times');
xlabel('$N_{\ell} ({\rm mm}^{-2})$','Interpreter','latex');
xticks([0 log10(3) 1]);
xticklabels({'1';'3';'$10^1$';'$10^2$'});
ylabel('pdf','Interpreter','latex');
yticks([])

n_l_H = random(n_l_0.base,[1 n]);
n_l_A = random(n_l_0.atherofact,[1 n]).*n_l_H;
h = histogram(log10(n_l_H),nBins,'FaceColor',col.Hbase,'BinLimits',[min(log10(n_l_H(:))) max(log10(n_l_H(:)))],'Normalization','pdf','Visible','off');
plot((h.BinEdges(1:end-1)+h.BinEdges(2:end))/2,smooth(h.Values,5),'linew',2,'color',col.Hbase);
h = histogram(log10(n_l_A),nBins,'FaceColor',col.Abase,'BinLimits',[min(log10(n_l_A(:))) max(log10(n_l_A(:)))],'Normalization','pdf','Visible','off');
plot((h.BinEdges(1:end-1)+h.BinEdges(2:end))/2,smooth(h.Values,5),'linew',2,'color',col.Abase,'marker','*','markerindices',5:5:nBins);

n_l_L_H = random(n_l.base,[1 n]);
n_l_L_A = random(n_l_0.atherofact,[1 n]).*n_l_L_H;
h = histogram(log10(n_l_L_H),nBins,'FaceColor',col.Hbase,'BinLimits',[min(log10(n_l_L_H(:))) max(log10(n_l_L_H(:)))],'Normalization','pdf','Visible','off');
plot((h.BinEdges(1:end-1)+h.BinEdges(2:end))/2,smooth(h.Values,5),'linew',2,'color',col.Hbase,'lines',ls.L);
h = histogram(log10(n_l_L_A),nBins,'FaceColor',col.Abase,'BinLimits',[min(log10(n_l_L_A(:))) max(log10(n_l_L_A(:)))],'Normalization','pdf','Visible','off');
plot((h.BinEdges(1:end-1)+h.BinEdges(2:end))/2,smooth(h.Values,5),'linew',2,'color',col.Abase,'lines',ls.L,'marker','*','markerindices',5:5:nBins);



%% lymphatic drainage rate
ax6 = nexttile([5 3]);hold all;
set(gca,'fontsize',15,'TickLabelInterpreter','latex');
text(-0.1,1,labels{4},'Units','normalized','FontSize',14,'fontweight','bold','fontname','times');
xlabel('$\tilde{q}_{\ell} ({\rm mm}^2/{\rm s})$','Interpreter','latex');
xlim(log10([1e-12 3e-11]))
xticks([-12 log10(3e-12) -11 log10(3e-11)]);
xticklabels({'$10^{-6}$';'$3\,10^{-6}$';'$10^{-5}$';'$3\,10^{-5}$'});
ylabel('pdf','Interpreter','latex');
yticks([])

q_l_H = 10.^random(q_l_0.base,[1 n]);
q_l_A = random(q_l_0.atherofact,[1 n]).*q_l_H;
line([log10(min(q_l_H)) log10(max(q_l_H))],1/(log10(max(q_l_H))-log10(min(q_l_H)))*[1 1],'color',col.Hbase,'linew',2);
h = histogram(log10(q_l_A),nBins,'FaceColor',col.Abase,'BinLimits',[min(log10(q_l_A(:))) max(log10(q_l_A(:)))],'Normalization','pdf','Visible','off');
plot((h.BinEdges(1:end-1)+h.BinEdges(2:end))/2,smooth(h.Values,5),'linew',2,'color',col.Abase,'marker','*','markerindices',5:5:nBins);

q_l_L_H = 10.^random(q_l.base,[1 n]);
q_l_L_A = random(q_l_0.atherofact,[1 n]).*q_l_L_H;
h = histogram(log10(q_l_L_H),nBins,'FaceColor',col.Hbase,'BinLimits',[min(log10(q_l_L_H(:))) max(log10(q_l_L_H(:)))],'Normalization','pdf','Visible','off');
plot((h.BinEdges(1:end-1)+h.BinEdges(2:end))/2,smooth(h.Values,5),'linew',2,'color',col.Hbase,'lines',ls.L);
h = histogram(log10(q_l_L_A),nBins,'FaceColor',col.Abase,'BinLimits',[min(log10(q_l_L_A(:))) max(log10(q_l_L_A(:)))],'Normalization','pdf','Visible','off');
plot((h.BinEdges(1:end-1)+h.BinEdges(2:end))/2,smooth(h.Values,5),'linew',2,'color',col.Abase,'lines',ls.L,'marker','*','markerindices',5:5:nBins);


if EXPORT
    exportgraphics(gcf,[homepath 'figures\distributions_post_Lambda.png'],'Resolution',300);
end