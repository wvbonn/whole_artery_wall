% radial model of interstitial fluid flow and solute transport in an
% endothelium-to-PVAT arterial wall with distributed vasa vasorum and
% lymphatic vessels 
% author: Willy Bonneuil, Université de Rennes 
% e-mail: willy.bonneuil@univ-rennes.fr

clear
close all
tic

currdir = pwd;
k = strfind(currdir,'\');
homepath = currdir(1:k(end));

ATHERO = 1==0;
N = 10000; 

savename = 'output_healthy_fluid_unfiltered_parameters.mat';

param = declare_fixed_parameters(ATHERO);

% parameters not affected by p_ap filtration + athero. factors
load('unfiltered_distributions.mat');
k_i_0 = k_i; k_a_0 = k_a; l_pv_0 = l_pv; n_v_0 = n_v; n_l_0 = n_l; q_l_0 = q_l;

for i = 1:N
    rng('shuffle');
    
    % varied parameters
    varpar(i).k_i = random(k_i.base,1);
    varpar(i).k_a = 10^random(k_a.base,1);
    varpar(i).t_p = random(t_p.base,1);
    varpar(i).p_lu = random(p_lu.base,1);
    varpar(i).l_pv = random(l_pv.base,1);
    varpar(i).n_v = random(n_v.base,1);
    varpar(i).n_l = random(n_l.base,1);
    varpar(i).q_l = 10^random(q_l.base,1);
    
    % simulations
    varpar(i).r = get_solution_grid(varpar(i),param);
    [p{i},u{i},Q_v{i},Q_l{i}] = radial_model_fluid(varpar(i),param,ATHERO);
    r{i} = varpar(i).r;
    [p_ia(i),p_ap(i)] = get_interface_pressures(p{i},r{i},param);
    u_i(i) = get_inner_layer_velocity(u{i},r{i},param);
    Q_v_l(i) = get_vv_fraction_lymph(Q_v{i},Q_l{i},r{i});
    [~,d_dil(i)] = get_dilution_distance(u{i},Q_v{i},r{i},param);
    
    if mod(i,100)==0; toc; end
end

save(savename,'r','p','u','Q_v','Q_l','u_i','p_ia','p_ap','Q_v_l','d_dil','param','varpar');

