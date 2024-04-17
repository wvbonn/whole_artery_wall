function q_lv = get_qlv_from_Lambda(v,p)

    % equivalent inner-layer/adventitia thickness
    tt = (p.t_im*p.t_a)/(p.t_im+p.t_a)*1e-3;
    % equivalent inner-layer/adventitia permeability    
    kk = v.k_im.*v.k_a./(v.k_im+v.k_a);
    % flux through lumen
    IM = kk/(p.mu/133).*v.p_lu./(2*tt.*p.r_eel*1e-3);
    % flux across vasa vasorum
    VV = v.l_vv*10*v.n_vv*p.l_ex*pi*p.d_vv*(p.p_vv-p.sigma*(p.pi_vv-p.pi_ref));
    % flux into lymphatics
    LV = v.Lambda*(VV+IM);
    q_lv = LV/v.n_lv/1e6;

end