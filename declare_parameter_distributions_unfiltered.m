clear
close all

savename = 'unfiltered_distributions.mat';

% inner layer permeability [m^2]
x_k_i = [0.2e-18 1.8e-18 100]; mu_k_i = mean(x_k_i(1:2)); s_k_i = (x_k_i(2)-x_k_i(1))/4;
k_i.base = trunc_norm_distr(mu_k_i,s_k_i,x_k_i);
x_k_i_a = [1 3.5 100]; mu_k_i_a = mean(x_k_i_a(1:2)); s_k_i_a = (x_k_i_a(2)-x_k_i_a(1))/4;
k_i.atherofact = trunc_norm_distr(mu_k_i_a,s_k_i_a,x_k_i_a);

% adventitial permeability [m^2]
x_k_a = [5e-18 5e-16 100]; 
k_a.base = trunc_unif_distr([log10(x_k_a(1:2)) 100]);
x_k_a_a = [2e-18 2e-16 100]; 
k_a.athero = trunc_unif_distr([log10(x_k_a_a(1:2)) 100]);

% PVAT thickness [mm]
x_t_p = [0.88 1.76 100];
t_p.base = trunc_unif_distr(x_t_p);

% lumen pressure [mmHg]
x_p_lu = [70 120 100];
p_lu.base = trunc_unif_distr(x_p_lu);

% vasa vasorum density [mm^-2]
x_n_v = [10 40 100]; mu_n_v = mean(x_n_v(1:2)); s_n_v = (x_n_v(2)-x_n_v(1))/4;
n_v.base = trunc_norm_distr(mu_n_v,s_n_v,x_n_v);
x_n_v_a = [2.05 3.05 100]; mu_n_v_a = mean(x_n_v_a(1:2)); s_n_v_a = (mu_n_v_a-x_n_v_a(1))/4;
n_v.atherofact = trunc_norm_distr(mu_n_v_a,s_n_v_a,x_n_v_a);

% vasa vasorum conductivity [cm/s/mmHg]
x_l_pv = [0.6e-7 1.4e-7 100]; mu_l_pv = mean(x_l_pv(1:2)); s_l_pv = (x_l_pv(2)-x_l_pv(1))/4;
l_pv.base = trunc_norm_distr(mu_l_pv,s_l_pv,x_l_pv);
x_l_pv_a = [1 3 100]; mu_l_pv_a = mean(x_l_pv_a(1:2)); s_l_pv_a = (mu_l_pv_a-x_l_pv_a(1))/4;
l_pv.atherofact = trunc_norm_distr(mu_l_pv_a,s_l_pv_a,x_l_pv_a);

% lymphatic density [mm^-2]
x_n_l = [0.5 5.5 100]; mu_n_l = mean(x_n_l(1:2)); s_n_l = (x_n_l(2)-x_n_l(1))/4;
n_l.base = trunc_norm_distr(mu_n_l,s_n_l,x_n_l);
x_n_l_a = [2 15 100]; nu_n_l_a = 7;
n_l.atherofact = trunc_chi2_distr(nu_n_l_a,x_n_l);

% lymphatic drainage [m^2/s]
x_q_l = [7e-12 2.8e-11 100];
q_l.base = trunc_unif_distr([log10(x_q_l(1:2)) 100]);
x_q_l_a = [1/8 1 100]; ndf_q_l_a = 8; ddf_q_l_a = 1; ct_q_l_a = 0.5;
q_l.atherofact = trunc_f_distr(x_q_l_a,ndf_q_l_a,ddf_q_l_a,ct_q_l_a);


save(savename,'k_i','k_a','t_p','p_lu','l_pv','n_l','n_v','q_l');


