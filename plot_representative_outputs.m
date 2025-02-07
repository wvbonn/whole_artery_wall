
clear
close all

currdir = pwd;
k = strfind(currdir,'\');
homepath = currdir(1:k(end));

EXPORT= 1==1;

param = declare_fixed_parameters(1==0);
r_ia = param.r_ia;
r_ap = r_ia+param.t_a;
D = 1e-5; % diffusivity [mm^2/s]

varpar.k_i = 1e-18;
varpar.k_a = 1e-17;
varpar.t_p = 1.3;
varpar.p_lu = 80;
varpar.l_pv = 1e-7;
varpar.n_v = 30;
varpar.n_l = 3;
varpar.q_l = 1e-11;
varpar.R_d = 30;

% simulation fluid
varpar.r = get_solution_grid(varpar,param);
[p,v,Q_v,Q_l] = radial_model_fluid(varpar,param,1==0);
r = varpar.r;
[p_ia,p_ap] = get_interface_pressures(p,r,param);
v_im = get_inner_layer_velocity(v,r,param);
Q_v_l = get_vv_fraction_lymph(Q_v,Q_l,r);
[dil,d_dil10] = get_dilution_distance(v,Q_v,r,param);

% normalise radius
r_norm = (r-min(r))/(max(r)-min(r));
% find indices for EEL and APV
id_ia = find(r>=r_ia,1);
id_ap = find(r>=r_ap,1);
r_ia_norm = r_norm(id_ia);
r_ap_norm = r_norm(id_ap);
r_p = r_norm(id_ia:end);
% normalise pressure
p_norm = p/(max(p)-min(p));
% normalise interstitial velocity
v_norm = (v-min(v))/(max(v)-min(v));
% normalise microvascular fluxes
Q_v_norm = Q_v(id_ia:end)/max(Q_l);
Q_l_norm = Q_l(id_ia:end)/max(Q_l);
% normalise dilution distance
d_dil10_norm = (d_dil10-min(r))/(max(r)-min(r));

figure('position',[50 50 1000 800],'color','w');
tiledlayout(3,2);
labels = get_subplot_labels('a':'z',6);

nexttile;hold on;
xlim([0 1])
ylim([0 1])
set(gca,'FontSize',14,'TickLabelInterpreter','latex');
ylabel({'velocity *'},'Interpreter','latex')
yticks([])
%xlabel('r','Interpreter','latex')
xticks([r_ia_norm r_ap_norm])
xticklabels({'$r_{ia}$';'$r_{ap}$'})
text(-0.15,1,labels{1},'Units','normalized','FontSize',13,'fontweight','bold','FontName','times');
plot(r_norm,v_norm,'linew',1.5,'color',0.5*[1 1 1]);
fill([r_norm(1) r_norm(id_ia) r_norm(id_ia) r_norm(1)],[v_norm(id_ia) v_norm(id_ia) v_norm(1) v_norm(1)],'k','facealpha',0,'edgecolor','k','linew',1.5);
text(r_norm(id_ia)+0.02,0.5*(v_norm(1)+v_norm(id_ia))+0.01,'average: $\bar{u}_{\rm im}$','interpreter','latex','fontsize',15);

nexttile;hold on;
xlim([0 1])
ylim([-0.1 1])
set(gca,'FontSize',14,'TickLabelInterpreter','latex');
ylabel({'pressure *'},'Interpreter','latex')
yticks([])
%xlabel('r','Interpreter','latex')
xticks([r_ia_norm r_ap_norm])
xticklabels({'$r_{ia}$';'$r_{ap}$'})
text(-0.15,1,labels{2},'Units','normalized','FontSize',13,'fontweight','bold','FontName','times');
plot(r_norm,p_norm,'linew',1.5,'color',0.5*[1 1 1]);
line([r_norm(1) r_norm(id_ia)],[p_norm(id_ia) p_norm(id_ia)],'lines','--','color','k');
scatter(r_norm(id_ia),p_norm(id_ia),30,'k','filled');
text(r_norm(id_ia),p_norm(id_ia)+0.1,'$p_{ia}$','Interpreter','latex','FontSize',15);
line([r_norm(1) r_norm(id_ap)],[p_norm(id_ap) p_norm(id_ap)],'lines','--','color','k');
scatter(r_norm(id_ap),p_norm(id_ap),30,'k','filled');
text(r_norm(id_ap),p_norm(id_ap)+0.1,'$p_{ap}$','Interpreter','latex','FontSize',15);

nexttile;hold on;
xlim([0 1.1])
ylim([min(Q_v_norm) 1])
set(gca,'FontSize',14,'TickLabelInterpreter','latex');
ylabel({'microvascular fluxes *'},'Interpreter','latex')
yticks([])
%xlabel('r','Interpreter','latex')
xticks([r_ia_norm r_ap_norm])
xticklabels({'$r_{ia}$';'$r_{ap}$'})
text(-0.15,1,labels{3},'Units','normalized','FontSize',13,'fontweight','bold','FontName','times');
fill([r_p fliplr(r_p)],[min(Q_v_norm)*ones(1,numel(r_p)) fliplr(Q_v_norm)],[1 0 0],'FaceAlpha',0.25,'linew',1.5,'edgecolor',[1 0 0]);
fill([r_p fliplr(r_p)],[min(Q_v_norm)*ones(1,numel(r_p)) fliplr(Q_l_norm)],[0 1 0],'FaceAlpha',0.25,'linew',1.5,'edgecolor',[0 1 0]);
plot(r_p,Q_v_norm,'linew',1.5,'color',0.5*[1 1 1]);
plot(r_p,Q_l_norm,'linew',1.5,'color',0.5*[1 1 1]);
text(0.9,1.05,'$q_{\ell}$','Interpreter','latex','FontSize',15);
text(0.05,0.94,'$Q_{\ell}=\int q_{\ell}(r) r {\rm d}r$','Interpreter','latex','FontSize',15);
text(0.9,0.95,'$q_v$','Interpreter','latex','FontSize',15);
text(0.4,0.5,'$Q_v=\int q_v(r) r {\rm d}r$','Interpreter','latex','FontSize',15);

nexttile;hold on;
xlim([0 1])
set(gca,'FontSize',14,'TickLabelInterpreter','latex');
ylabel({'proportion of';'luminal fluid *'},'Interpreter','latex')
%xlabel('r','Interpreter','latex')
xticks([r_ia_norm r_ap_norm])
xticklabels({'$r_{ia}$';'$r_{ap}$'})
text(-0.15,1,labels{4},'Units','normalized','FontSize',13,'fontweight','bold','FontName','times');
plot(r_p,dil(id_ia:end),'linew',1.5,'color',0.5*[1 1 1]);
line([r_norm(1) d_dil10_norm],[0.1 0.1],'lines','--','color','k');
line(d_dil10_norm*[1 1],[0 0.1],'lines','--','color','k');
scatter(d_dil10_norm,0.1,30,'black','filled');
text(d_dil10_norm+0.01,0.3,{'10\%';'dilution';'distance'},'Interpreter','latex','FontSize',14);


%% transport outputs
% simulation 
id = id_ia:numel(r_norm);
varpar.r = varpar.r(id);
param.D = D;
[c,dc] = radial_model_tsp(1e-3*v(id),Q_v(id),Q_l(id),varpar,param);
if any(((dc>-0.004)+(c<1))>1) % concentration gradient low and concentration low so as not to count the cases where there are zero-gradient points because of advection-driven accumulation
    t_dc_ccl = (varpar.r(min(intersect(find(dc>-0.004),find(c<1))))-(r_ia+param.t_a))/varpar.t_p;
else
    t_dc_ccl = 1;
end
r = varpar.r;

% normalise radius
r_norm = (r-min(r))/(max(r)-min(r));
r_ia_norm = 0;
r_ap_norm = (r_ia+param.t_a-min(r))/(max(r)-min(r));
% find indices for t_dc_ccl. 
id_tdc = min(intersect(find(dc>-0.004),find(c<1)));
% normalise concentration gradient
dc_norm = (dc-max(dc))/(max(dc)-min(dc));

nexttile;hold on;
xlim([0 1])
ylim([0 1])
set(gca,'FontSize',14,'TickLabelInterpreter','latex');
ylabel({'solute';'concentration *'},'Interpreter','latex')
%xlabel('r','Interpreter','latex')
xticks([r_ia_norm r_ap_norm])
xticklabels({'$r_{ia}$';'$r_{ap}$'})
text(-0.15,1,labels{5},'Units','normalized','FontSize',13,'fontweight','bold','FontName','times');
plot(r_norm,c,'linew',1.5,'color',0.5*[1 1 1]);

nexttile;hold on;
xlim([0 1])
set(gca,'FontSize',14,'TickLabelInterpreter','latex');
ylabel({'solute';'gradient *'},'Interpreter','latex')
yticks(0)
%xlabel('r','Interpreter','latex')
xticks([r_ia_norm r_ap_norm])
xticklabels({'$r_{ia}$';'$r_{ap}$'})
text(-0.15,1,labels{6},'Units','normalized','FontSize',13,'fontweight','bold','FontName','times');
plot(r_norm,dc_norm,'linew',1.5,'color',0.5*[1 1 1]);
text(0.01,-0.07,'-0.4\%','Interpreter','latex','FontSize',15);
line([0 r_norm(id_tdc)],[-0.004 -0.004],'lines','--','color','k');
line([r_norm(id_tdc) r_norm(id_tdc)],[min(dc_norm) -0.004],'lines','--','color','k');
scatter(r_norm(id_tdc),-0.004,30,"black","filled");
annotation('textbox',[.83 .06 .1 .1],'string','$t^*_{\rm DC,PVAT}$','Interpreter','latex','FontSize',15,'edgecolor','none','VerticalAlignment','bottom');
yl = get(gca,'YLim');
scatter(0,yl(1),30,'k','filled');
annotation('textbox',[.48 .09 .1 .1],'string','${\rm max} ||\nabla C||$','fontsize',15,'interpreter','latex','edgecolor','none','VerticalAlignment','bottom');

if EXPORT
    exportgraphics(gcf,[homepath 'figures\representative_model_outputs.png'],'Resolution',300);
end