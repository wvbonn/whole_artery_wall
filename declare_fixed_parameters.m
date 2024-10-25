function param = declare_fixed_parameters(ATHERO)

    param.rho = 1050; % density of interstitial fluid (taken like lymph)
    param.mu = 1.1e-3; % dynamic viscosity of water at 37°C (SI)
    param.r_ia = 2; % EEL radius [mm]
    if ATHERO
        param.t_i = 0.78; % intimal-medial-adventitial thickness [mm]
        param.t_a = 0.92; % adventitia thickness [mm]
    else
        param.t_i = 0.34; % intimal-medial-adventitial thickness [mm]
        param.t_a = 0.56; % adventitia thickness [mm]
    end
    param.k_p = 1e-15; % PVAT permeability [m^2]
    param.sigma = 0.8; % osmotic reflection coefficient, assumed like lymph node
    param.p_v = 25; % hydrostatic pressure, vasa vasorum [mmHg]
    param.pi_v = 20; % osmotic pressure, vasa vasorum [mmHg]
    param.pi_ref = 8; % osmotic pressure, tissue [mmHg]
    param.d_v = 10e-3; % baseline capillary diameter [mm]
    param.c_0 = 1; % inlet concentration [-]

    param.s = 0.01; % mesh size [mm]
    
end