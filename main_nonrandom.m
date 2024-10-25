% radial model of interstitial fluid flow and solute transport in an
% endothelium-to-PVAT arterial wall with distributed vasa vasorum and
% lymphatic vessels 
% author: Willy Bonneuil, Université de Rennes 
% e-mail: willy.bonneuil@univ-rennes.fr

clear
close all

ATHERO = 1==0;
TSP = 1==1; % solve for solute transport (true) or just for fluid flow (false)

if ATHERO
    if TSP
        savename = 'nonrnd_athero_tsp.mat';
    end
else
    if TSP
        savename = 'nonrnd_healthy_tsp.mat';
    else
        savename = 'nonrnd_healthy_fluid.mat';
    end
end

param = declare_fixed_parameters(ATHERO);
D = 1e6*[1e-10 1e-11 1e-12]; % diffusivities [mm^2/s]

t_pv = 0.88*ones(1,8);
p_lu = 90*ones(1,8);
if ~TSP
    k_i = repmat([0.7 1.2],1,4)*1e-18; 
    k_a = repmat([1 1 10 10],1,2)*1e-17; 
    l_pv = 7.5e-8*ones(1,8); 
    n_v = 25*ones(1,8); 
    n_l = 3*ones(1,8);
    q_l_0 = [0.8*ones(1,4) 0.9*ones(1,4)]*1e-11;
else
    k_a = 5e-17*ones(1,8);
    R_d = repmat([3 3 30 30],1,2); 
    if ATHERO
        k_i = repmat([1.2 2.8],1,4)*1e-18; 
        l_pv = 2*7.5e-8*ones(1,8); 
        n_v = 60*ones(1,8); 
        n_l = 10*ones(1,8);
        q_l_0 = [0.2*ones(1,4) 1*ones(1,4)]*1e-11;
    else
        k_i = repmat([0.7 1.2],1,4)*1e-18; 
        l_pv = 7.5e-8*ones(1,8); 
        n_v = 25*ones(1,8); 
        n_l = 3*ones(1,8);
        q_l_0 = [0.8*ones(1,4) 0.9*ones(1,4)]*1e-11;
    end
end

for i = 1:size(k_i,2)

    varpar(i).k_i = k_i(i);
    varpar(i).k_a = k_a(i);
    varpar(i).t_p = t_pv(i);
    varpar(i).p_lu = p_lu(i);
    varpar(i).l_pv = l_pv(i);
    varpar(i).n_v = n_v(i);
    varpar(i).n_l = n_l(i);
    varpar(i).q_l = q_l_0(i);
    if TSP
        varpar(i).R_d = R_d(i);
    end

    % simulations fluid
    varpar(i).r = get_solution_grid(varpar(i),param);
    [p{i},u{i},Q_v{i},Q_l{i}] = radial_model_fluid(varpar(i),param,ATHERO);
    if ~TSP
        r{i} = varpar(i).r;
    else
        r{i} = varpar(i).r(varpar(i).r>=param.r_ia);
    end

    % simulations transport
    if TSP
        id = find(varpar(i).r>=param.r_ia);
        varpar(i).r = varpar(i).r(id);
        for j = 1:numel(D)
            Pe(i,j) = varpar(i).k_i*1e6/(param.mu/133)*varpar(i).p_lu/param.t_i*(param.t_a+varpar(i).t_p)/D(j);
            Da_v(i,j) = varpar(i).l_pv*10*varpar(i).n_v*param.d_v*pi*(param.p_v-param.sigma*(param.pi_v-param.pi_ref))*(param.t_a+varpar(i).t_p)^2/D(j);
            Da_l(i,j) = varpar(i).n_l*varpar(i).q_l*1e6*(param.t_a+varpar(i).t_p)^2/D(j);
            param.D = D(j);
            [c{i,j},dc{i,j}] = radial_model_tsp(1e-3*u{i}(id),Q_v{i}(id),Q_l{i}(id),varpar(i),param);
        end
    end

end


if TSP
    save(savename,'r','Pe','Da_v','Da_l','c','dc','param','varpar');
else
    save(savename,'r','p','u','Q_v','Q_l','param','varpar');
end

