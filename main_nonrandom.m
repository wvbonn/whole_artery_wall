% radial model of interstitial flow and solute transport

clear
close all

homepath = 'C:\Users\Willy\Work\PhD\ATLO\clean\';
ATHERO = 1==1;
TSP = 1==1; % solve for solute transport (true) or just for fluid flow (false)

if ATHERO
    if TSP
        savename = 'nonrnd_1D_athero_tsp.mat';
    end
else
    if TSP
        savename = 'nonrnd_1D_healthy_tsp.mat';
    else
        savename = 'nonrnd_1D_healthy_fluid.mat';
    end
end

param = declare_fixed_parameters_1D(ATHERO);
D = 1e6*[1e-10 1e-11 1e-12]; % diffusivities [mm^2/s]

t_pv = 0.88*ones(1,8);
p_lu = 90*ones(1,8);
if ~TSP
    k_im = repmat([0.6 1.5],1,4)*1e-18; 
    k_a = repmat([1 1 10 10],1,2)*1e-17; 
    l_vv = 7.5e-8*ones(1,8); 
    n_vv = 25*ones(1,8); 
    n_lv = 3*ones(1,8);
    q_lv_0 = [0.9*ones(1,4) 1.2*ones(1,4)]*1e-11;
else
    k_a = 5e-17*ones(1,8);
    R_d = repmat([3 3 30 30],1,2); 
    if ATHERO
        k_im = repmat([1.2 2.8],1,4)*1e-18; 
        l_vv = 2*7.5e-8*ones(1,8); 
        n_vv = 60*ones(1,8); 
        n_lv = 10*ones(1,8);
        q_lv_0 = [0.2*ones(1,4) 1*ones(1,4)]*1e-11;
    else
        k_im = repmat([0.6 1.5],1,4)*1e-18; 
        l_vv = 7.5e-8*ones(1,8); 
        n_vv = 25*ones(1,8); 
        n_lv = 3*ones(1,8);
        q_lv_0 = [0.9*ones(1,4) 1.2*ones(1,4)]*1e-11;
    end
end

for i = 1:size(k_im,2)

    varpar(i).k_im = k_im(i);
    varpar(i).k_a = k_a(i);
    varpar(i).t_pv = t_pv(i);
    varpar(i).p_lu = p_lu(i);
    varpar(i).l_vv = l_vv(i);
    varpar(i).n_vv = n_vv(i);
    varpar(i).n_lv = n_lv(i);
    varpar(i).q_lv = q_lv_0(i);
    varpar(i).R_d = R_d(i);

    % simulations fluid
    varpar(i).r = get_solution_grid(varpar(i),param);
    [p{i},v{i},Q_vv{i},Q_lv{i}] = radial_model_fluid(varpar(i),param,ATHERO);
    if ~TSP
        r{i} = varpar(i).r;
    else
        r{i} = varpar(i).r(varpar(i).r>=param.r_eel);
    end

    % simulations transport
    if TSP
        id = find(varpar(i).r>=param.r_eel);
        varpar(i).r = varpar(i).r(id);
        for j = 1:numel(D)
            Pe(i,j) = varpar(i).k_im*1e6/(param.mu/133)*varpar(i).p_lu/param.t_im*(param.t_a+varpar(i).t_pv)/D(j);
            Da_vv(i,j) = varpar(i).l_vv*10*varpar(i).n_vv*param.d_vv*param.l_ex*pi*(param.p_vv-param.sigma*(param.pi_vv-param.pi_ref))*(param.t_a+varpar(i).t_pv)^2/D(j);
            Da_lv(i,j) = varpar(i).n_lv*varpar(i).q_lv*1e6*(param.t_a+varpar(i).t_pv)^2/D(j);
            param.D = D(j);
            [c{i,j},dc{i,j}] = radial_model_tsp(1e-3*v{i}(id),Q_vv{i}(id),Q_lv{i}(id),varpar(i),param);
        end
    end

end


if TSP
    save(savename,'r','Pe','Da_vv','Da_lv','c','dc','param','varpar');
else
    save(savename,'r','p','v','Q_vv','Q_lv','param','varpar');
end

