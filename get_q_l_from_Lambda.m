function q_l = get_q_l_from_Lambda(v,p)

    % radii [mm]
    r_lu = p.r_ia-p.t_i;
    r_ap = p.r_ia+p.t_a;
    % equivalent inner-layer/adventitia permeability    
    kk = log(r_ap/r_lu)/(log(p.r_ia/r_lu)/v.k_i+log(r_ap/p.r_ia)/v.k_a);
    % mass through lumen
    I = kk/(p.mu/133).*(v.p_lu+2.5)*2./(v.t_p^2+2*r_ap*v.t_p)*1e6;
    % mass across vasa vasorum
    V = v.l_pv*10*v.n_v*pi*p.d_v*(p.p_v+2.5-p.sigma*(p.pi_v-p.pi_ref));
    % mass into lymphatics
    L = v.Lambda*(V+I);
    q_l = L/v.n_l/1e6;


end