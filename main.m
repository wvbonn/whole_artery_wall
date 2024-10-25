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
TSP = 1==1; % solve for solute transport (true) or just for fluid flow (false)
N = 1000; 

if ATHERO
    if TSP
        savename = 'output_athero_tsp.mat';
    else
        savename = 'output_athero_fluid.mat';
    end
else
    if TSP
        savename = 'output_healthy_tsp.mat';
    else
        savename = 'output_healthy_fluid.mat';
    end
end

param = declare_fixed_parameters(ATHERO);
D = 1e6*[1e-10 1e-11 1e-12]; % diffusivities [mm^2/s]
L = [0.9 1.6]; % bounds of physiological p_ap filtration

% parameters not affected by p_apv filtration + athero. factors
load('unfiltered_distributions.mat');
k_i_0 = k_i; k_a_0 = k_a; l_pv_0 = l_pv; n_v_0 = n_v; n_l_0 = n_l; q_l_0 = q_l;
% parameters affected by p_apv filtration
load(['filtered_Lambda_' num2str(L(1)) '_' num2str(L(2)) '.mat']);

for i = 1:N
    rng('shuffle');
    
    % varied parameters
    varpar(i).k_i = random(k_i.base,1);
    if ATHERO; varpar(i).k_i = varpar(i).k_i*random(k_i_0.atherofact,1); end
    if ATHERO; varpar(i).k_a = 10^random(k_a_0.athero,1); else; varpar(i).k_a = random(k_a.base,1); end
    varpar(i).t_p = random(t_p.base,1);
    varpar(i).p_lu = random(p_lu.base,1);
    varpar(i).l_pv = random(l_pv.base,1);
    if ATHERO; varpar(i).l_pv = varpar(i).l_pv*random(l_pv_0.atherofact,1); end
    varpar(i).n_v = random(n_v.base,1);
    if ATHERO; varpar(i).n_v = varpar(i).n_v*random(n_v_0.atherofact,1); end
    varpar(i).n_l = random(n_l.base,1);
    if ATHERO; varpar(i).n_l = varpar(i).n_l*random(n_l_0.atherofact,1); end
    varpar(i).Lambda = unifrnd(L(1),L(2),1);
    varpar(i).q_l = get_q_l_from_Lambda(varpar(i),param);
    if ATHERO; varpar(i).q_l = varpar(i).q_l*random(q_l_0.atherofact,1); end
    if TSP; varpar(i).R_d = 10^unifrnd(0,2,1); end
 
    % simulations fluid
    varpar(i).r = get_solution_grid(varpar(i),param);
    [p{i},u{i},Q_v{i},Q_l{i}] = radial_model_fluid(varpar(i),param,ATHERO);
    r{i} = varpar(i).r;
    if ~TSP
        [p_ia(i),p_ap(i)] = get_interface_pressures(p{i},r{i},param);
        u_i(i) = get_inner_layer_velocity(u{i},r{i},param);
        Q_v_l(i) = get_vv_fraction_lymph(Q_v{i},Q_l{i},r{i});
        [~,d_dil(i)] = get_dilution_distance(u{i},Q_v{i},r{i},param);
    end

    % simulations transport. restrict pressure & velocity fields to
    % adventitia + PVAT
    if TSP
        id = find(varpar(i).r>=param.r_ia);
        r{i} = varpar(i).r(id);
        varpar(i).r = varpar(i).r(id);
        for j = 1:numel(D)
            Pe(i,j) = varpar(i).k_i*1e6/(param.mu/133)*varpar(i).p_lu/param.t_i*(param.t_a+varpar(i).t_p)/D(j);
            Da_v(i,j) = varpar(i).l_pv*10*varpar(i).n_v*param.d_v*pi*(param.p_v-param.sigma*(param.pi_v-param.pi_ref))*(param.t_a+varpar(i).t_p)^2/D(j);
            Da_l(i,j) = varpar(i).n_l*varpar(i).q_l*1e6*(param.t_a+varpar(i).t_p)^2/D(j);
            param.D = D(j);
            [c{i,j},dc{i,j}] = radial_model_tsp(1e-3*u{i}(id),Q_v{i}(id),Q_l{i}(id),varpar(i),param);
            % max gradient and DC-CCL19 transport distance
            [max_dc(i,j),max_dc_id(i,j)] = max(abs(dc{i,j}));
            if any(((dc{i,j}>-0.004)+(c{i,j}<0.9))>1) % concentration gradient low and concentration low
                t_dc_ccl(i,j) = (r{i}(min(intersect(find(dc{i,j}>-0.004),find(c{i,j}<0.9))))-(param.r_ia+param.t_a))/varpar(i).t_p;
            else
                t_dc_ccl(i,j) = 1;
            end
        end
    end
    
    if mod(i,100)==0; toc; end
end



if TSP
    save(savename,'r','Pe','Da_v','Da_l','c','dc','max_dc','max_dc_id','t_dc_ccl','param','varpar');
else
    save(savename,'r','p','u','Q_v','Q_l','u_i','p_ia','p_ap','Q_v_l','d_dil','param','varpar');
end

