% schematics explaining the point outputs

clear
close all

homepath = 'C:\Users\Willy\Work\PhD\ATLO\clean\';
EXPORT= 1==1;

param = declare_fixed_parameters_1D(1==0);
r_eel = param.r_eel;
r_apv = r_eel+param.t_a;
D = 1e-5; % diffusivity [mm^2/s]

varpar.k_im = 1e-18;
varpar.k_a = 1e-17;
varpar.t_pv = 1.3;
varpar.p_lu = 80;
varpar.l_vv = 1e-7;
varpar.n_vv = 30;
varpar.n_lv = 3;
varpar.q_lv = 1e-11;
varpar.R_d = 30;

% simulation fluid
varpar.r = get_solution_grid(varpar,param);
[p,v,Q_vv,Q_lv] = radial_model_fluid(varpar,param,1==0);
r = varpar.r;
[p_ma,p_apv] = get_interface_pressures(p,r,param);
v_im = get_inner_layer_velocity(v,r,param);
Q_vv_lv = get_vv_fraction_lymph(Q_vv,Q_lv,r);
[dil,d_dil10] = get_dilution_distance(v,Q_vv,r,param);

% normalise radius
r_norm = (r-min(r))/(max(r)-min(r));
% find indices for EEL and APV
id_eel = find(r>=r_eel,1);
id_apv = find(r>=r_apv,1);
r_eel_norm = r_norm(id_eel);
r_apv_norm = r_norm(id_apv);
r_peri = r_norm(id_eel:end);
% normalise pressure
p_norm = p/(max(p)-min(p));
% normalise interstitial velocity
v_norm = (v-min(v))/(max(v)-min(v));
% normalise microvascular fluxes
Q_vv_norm = Q_vv(id_eel:end)/max(Q_lv);
Q_lv_norm = Q_lv(id_eel:end)/max(Q_lv);
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
xlabel('r','Interpreter','latex')
xticks([r_eel_norm r_apv_norm])
xticklabels({'m-a';'a-pv'})
text(-0.15,1,labels{1},'Units','normalized','FontSize',13,'fontweight','bold','FontName','times');
plot(r_norm,v_norm,'linew',1.5,'color',0.5*[1 1 1]);
fill([r_norm(1) r_norm(id_eel) r_norm(id_eel) r_norm(1)],[v_norm(id_eel) v_norm(id_eel) v_norm(1) v_norm(1)],'k','facealpha',0,'edgecolor','k','linew',1.5);
text(r_norm(id_eel)+0.02,0.5*(v_norm(1)+v_norm(id_eel))+0.01,'average: $\bar{v}_{\rm im}$','interpreter','latex','fontsize',15);

nexttile;hold on;
xlim([0 1])
ylim([-0.1 1])
set(gca,'FontSize',14,'TickLabelInterpreter','latex');
ylabel({'pressure *'},'Interpreter','latex')
yticks([])
xlabel('r','Interpreter','latex')
xticks([r_eel_norm r_apv_norm])
xticklabels({'m-a';'a-pv'})
text(-0.15,1,labels{2},'Units','normalized','FontSize',13,'fontweight','bold','FontName','times');
plot(r_norm,p_norm,'linew',1.5,'color',0.5*[1 1 1]);
line([r_norm(1) r_norm(id_eel)],[p_norm(id_eel) p_norm(id_eel)],'lines','--','color','k');
scatter(r_norm(id_eel),p_norm(id_eel),30,'k','filled');
text(r_norm(id_eel),p_norm(id_eel)+0.1,'$p_{\rm ma}$','Interpreter','latex','FontSize',15);
line([r_norm(1) r_norm(id_apv)],[p_norm(id_apv) p_norm(id_apv)],'lines','--','color','k');
scatter(r_norm(id_apv),p_norm(id_apv),30,'k','filled');
text(r_norm(id_apv),p_norm(id_apv)+0.1,'$p_{\rm apv}$','Interpreter','latex','FontSize',15);

nexttile;hold on;
xlim([0 1.1])
ylim([min(Q_vv_norm) 1])
set(gca,'FontSize',14,'TickLabelInterpreter','latex');
ylabel({'microvascular fluxes *'},'Interpreter','latex')
yticks([])
xlabel('r','Interpreter','latex')
xticks([r_eel_norm r_apv_norm])
xticklabels({'m-a';'a-pv'})
text(-0.15,1,labels{3},'Units','normalized','FontSize',13,'fontweight','bold','FontName','times');
fill([r_peri fliplr(r_peri)],[min(Q_vv_norm)*ones(1,numel(r_peri)) fliplr(Q_vv_norm)],[1 0 0],'FaceAlpha',0.25,'linew',1.5,'edgecolor',[1 0 0]);
fill([r_peri fliplr(r_peri)],[min(Q_vv_norm)*ones(1,numel(r_peri)) fliplr(Q_lv_norm)],[0 1 0],'FaceAlpha',0.25,'linew',1.5,'edgecolor',[0 1 0]);
plot(r_peri,Q_vv_norm,'linew',1.5,'color',0.5*[1 1 1]);
plot(r_peri,Q_lv_norm,'linew',1.5,'color',0.5*[1 1 1]);
text(0.9,1.05,'$q_{\rm lv}$','Interpreter','latex','FontSize',15);
text(0.05,0.94,'$Q_{\rm lv}=\int q_{\rm lv}(r) r {\rm d}r$','Interpreter','latex','FontSize',15);
text(0.9,0.95,'$q_{\rm vv}$','Interpreter','latex','FontSize',15);
text(0.4,0.5,'$Q_{\rm vv}=\int q_{\rm vv}(r) r {\rm d}r$','Interpreter','latex','FontSize',15);

nexttile;hold on;
xlim([0 1])
set(gca,'FontSize',14,'TickLabelInterpreter','latex');
ylabel({'proportion of';'luminal fluid *'},'Interpreter','latex')
xlabel('r','Interpreter','latex')
xticks([r_eel_norm r_apv_norm])
xticklabels({'m-a';'a-pv'})
text(-0.15,1,labels{4},'Units','normalized','FontSize',13,'fontweight','bold','FontName','times');
plot(r_peri,dil(id_eel:end),'linew',1.5,'color',0.5*[1 1 1]);
line([r_norm(1) d_dil10_norm],[0.1 0.1],'lines','--','color','k');
line(d_dil10_norm*[1 1],[0 0.1],'lines','--','color','k');
scatter(d_dil10_norm,0.1,30,'black','filled');
text(d_dil10_norm+0.01,0.3,{'10\%';'dilution';'distance'},'Interpreter','latex','FontSize',14);


%% transport outputs
% simulation 
id = id_eel:numel(r_norm);
varpar.r = varpar.r(id);
param.D = D;
[c,dc] = radial_model_tsp(1e-3*v(id),Q_vv(id),Q_lv(id),varpar,param);
if any(((dc>-0.004)+(c<1))>1) % concentration gradient low and concentration low so as not to count the cases where there are zero-gradient points because of advection-driven accumulation
    t_dc_ccl = (varpar.r(min(intersect(find(dc>-0.004),find(c<1))))-(param.r_eel+param.t_a))/varpar.t_pv;
else
    t_dc_ccl = 1;
end
r = varpar.r;

% normalise radius
r_norm = (r-min(r))/(max(r)-min(r));
r_eel_norm = 0;
r_apv_norm = (param.r_eel+param.t_a-min(r))/(max(r)-min(r));
% find indices for t_dc_ccl. 
id_tdc = min(intersect(find(dc>-0.004),find(c<1)));
% normalise concentration gradient
dc_norm = (dc-max(dc))/(max(dc)-min(dc));

nexttile;hold on;
xlim([0 1])
ylim([0 1])
set(gca,'FontSize',14,'TickLabelInterpreter','latex');
ylabel({'solute';'concentration *'},'Interpreter','latex')
xlabel('r','Interpreter','latex')
xticks([r_eel_norm r_apv_norm])
xticklabels({'m-a';'a-pv'})
text(-0.15,1,labels{5},'Units','normalized','FontSize',13,'fontweight','bold','FontName','times');
plot(r_norm,c,'linew',1.5,'color',0.5*[1 1 1]);

nexttile;hold on;
xlim([0 1])
set(gca,'FontSize',14,'TickLabelInterpreter','latex');
ylabel({'solute';'gradient *'},'Interpreter','latex')
yticks(0)
xlabel('r','Interpreter','latex')
xticks([r_eel_norm r_apv_norm])
xticklabels({'m-a';'a-pv'})
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