% p_ap against Lambda. healthy config, baseline
% distributions

clear
close all

currdir = pwd;
k = strfind(currdir,'\');
homepath = currdir(1:k(end));

EXPORT = 1==1;

fname = 'output_healthy_fluid_unfiltered_parameters.mat';
param = declare_fixed_parameters(1==0);
col.L = 0.4*[1 1 1];
ls.L = '--';

% load outputs
load(fname);
n = numel(p_ap);
L = get_Lambda(varpar,param,n);
save(fname,'L','-append');

% plot
figure('position',[100 100 600 350],'color','w');
hold all;
set(gca,'fontsize',15,'TickLabelInterpreter','latex');
xlabel('$\Lambda$','Interpreter','latex'); % mass-flux balance
ylabel('$p_{ap}$ (mmHg)','Interpreter','latex');
text(0.8,30,'healthy configuration, initial distributions','fontsize',14,'interpreter','latex');
xlim([0 5])
ylim([-90 30])
scatter(L,p_ap,1);
xl = get(gca,'xlim'); yl = get(gca,'ylim');
plot([xl(1) xl(2)],[-10 -10],'linew',1.5,'color',[0.8 0 0]);
plot([xl(1) xl(2)],[5 5],'linew',1.5,'color',[0.8 0 0]);
text(2,-3,'{\bf physiologic interval}','fontsize',14,'Interpreter','latex');
text(1.25,-88,'{\bf filtration bounds}','fontsize',14,'interpreter','latex','rotation',90);
plot(0.9*[1 1],[yl(1) yl(1)+0.95*(yl(2)-yl(1))],'linew',1.5,'color',col.L,'lines','--');
plot(1.6*[1 1],[yl(1) yl(1)+0.95*(yl(2)-yl(1))],'linew',1.5,'color',col.L,'lines','--');


if EXPORT
    exportgraphics(gcf,[homepath 'figures\Lambda_p_ap.png'],'Resolution',300);
end


function L = get_Lambda(varpar,param,n)
    % equivalent inner-layer/adventitia thickness
    r_ap = param.r_ia+param.t_a;
    for i = 1:n
        % equivalent inner-layer/adventitia permeability    
        kk = 1/(0.36/varpar(i).k_i+0.63/varpar(i).k_a);
        % flux through lumen
        IM = kk/(param.mu/133).*(varpar(i).p_lu+2.5)*2./(varpar(i).t_p^2+2*r_ap*varpar(i).t_p)*1e6;
        % flux across vasa vasorum
        VV = varpar(i).l_pv*10*varpar(i).n_v*pi*param.d_v*(param.p_v+2.5-param.sigma*(param.pi_v-param.pi_ref));
        % flux into lymphatics
        LV = varpar(i).n_l*1e6*varpar(i).q_l;
        L(i) = LV/(VV+IM);
    end
end




