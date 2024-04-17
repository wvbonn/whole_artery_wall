% radial model of interstitial flow and solute transport

clear
close all
tic

homepath = 'C:\Users\Willy\Work\PhD\ATLO\clean\';
ATHERO = 1==1;
TSP = 1==1; % solve for solute transport (true) or just for fluid flow (false)
N = 1000; 

if ATHERO
    if TSP
        savename = 'output_1D_athero_tsp.mat';
    else
        savename = 'output_1D_athero_fluid.mat';
    end
else
    if TSP
        savename = 'output_1D_healthy_tsp.mat';
    else
        savename = 'output_1D_healthy_fluid.mat';
    end
end

param = declare_fixed_parameters_1D(ATHERO);
D = 1e6*[1e-10 1e-11 1e-12]; % diffusivities [mm^2/s]
L = [0.7 1.3]; % bounds of physiological p_apv filtration

% parameters not affected by p_apv filtration + athero. factors
load([homepath 'distributions\baseline.mat']);
k_im_0 = k_im; k_a_0 = k_a; l_vv_0 = l_vv; n_vv_0 = n_vv; n_lv_0 = n_lv; q_lv_0 = q_lv;
% parameters affected by p_apv filtration
load([homepath 'distributions\filtered_L_' num2str(L(1)) '_' num2str(L(2)) '.mat']);

for i = 1:N
    rng('shuffle');
    
    % varied parameters
    varpar(i).k_im = random(k_im.base,1);
    if ATHERO; varpar(i).k_im = varpar(i).k_im*random(k_im_0.atherofact,1); end
    if ATHERO; varpar(i).k_a = 10^random(k_a_0.athero,1); else; varpar(i).k_a = random(k_a.base,1); end
    varpar(i).t_pv = random(t_pv.base,1);
    varpar(i).p_lu = random(p_lu.base,1);
    varpar(i).l_vv = random(l_vv.base,1);
    if ATHERO; varpar(i).l_vv = varpar(i).l_vv*random(l_vv_0.atherofact,1); end
    varpar(i).n_vv = random(n_vv.base,1);
    if ATHERO; varpar(i).n_vv = varpar(i).n_vv*random(n_vv_0.atherofact,1); end
    varpar(i).n_lv = random(n_lv.base,1);
    if ATHERO; varpar(i).n_lv = varpar(i).n_lv*random(n_lv_0.atherofact,1); end
    varpar(i).Lambda = unifrnd(L(1),L(2),1);
    varpar(i).q_lv = get_qlv_from_Lambda(varpar(i),param);
    if ATHERO; varpar(i).q_lv = varpar(i).q_lv*random(q_lv_0.atherofact,1); end
    if TSP; varpar(i).R_d = 10^unifrnd(0,2,1); end
 
    % simulations fluid
    varpar(i).r = get_solution_grid(varpar(i),param);
    [p{i},v{i},Q_vv{i},Q_lv{i}] = radial_model_fluid(varpar(i),param,ATHERO);
    r{i} = varpar(i).r;
    if ~TSP
        [p_ma(i),p_apv(i)] = get_interface_pressures(p{i},r{i},param);
        v_im(i) = get_inner_layer_velocity(v{i},r{i},param);
        Q_vv_lv(i) = get_vv_fraction_lymph(Q_vv{i},Q_lv{i},r{i});
        [~,d_dil(i)] = get_dilution_distance(v{i},Q_vv{i},r{i},param);
    end

    % simulations transport. restrict pressure & velocity fields to
    % adventitia + PVAT
    if TSP
        id = find(varpar(i).r>=param.r_eel);
        r{i} = varpar(i).r(id);
        varpar(i).r = varpar(i).r(id);
        for j = 1:numel(D)
            Pe(i,j) = varpar(i).k_im*1e6/(param.mu/133)*varpar(i).p_lu/param.t_im*(param.t_a+varpar(i).t_pv)/D(j);
            Da_vv(i,j) = varpar(i).l_vv*10*varpar(i).n_vv*param.d_vv*param.l_ex*pi*(param.p_vv-param.sigma*(param.pi_vv-param.pi_ref))*(param.t_a+varpar(i).t_pv)^2/D(j);
            Da_lv(i,j) = varpar(i).n_lv*varpar(i).q_lv*1e6*(param.t_a+varpar(i).t_pv)^2/D(j);
            param.D = D(j);
            [c{i,j},dc{i,j}] = radial_model_tsp(1e-3*v{i}(id),Q_vv{i}(id),Q_lv{i}(id),varpar(i),param);
            % max gradient and DC-CCL19 transport distance
            [max_dc(i,j),max_dc_id(i,j)] = max(abs(dc{i,j}));
            if any(((dc{i,j}>-0.004)+(c{i,j}<1))>1) % concentration gradient low and concentration low so as not to count the cases where there are zero-gradient points because of advection-driven accumulation
                t_dc_ccl(i,j) = (r{i}(min(intersect(find(dc{i,j}>-0.004),find(c{i,j}<1))))-(param.r_eel+param.t_a))/varpar(i).t_pv;
            else
                t_dc_ccl(i,j) = 1;
            end
        end
    end
    
    if mod(i,100)==0; toc; end
end

if TSP
    save(savename,'r','Pe','Da_vv','Da_lv','c','dc','max_dc','max_dc_id','t_dc_ccl','param','varpar');
else
    save(savename,'r','p','v','Q_vv','Q_lv','v_im','p_ma','p_apv','Q_vv_lv','d_dil','param','varpar');
end

