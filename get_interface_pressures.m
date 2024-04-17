function [p_ma,p_apv] = get_interface_pressures(p,r,param)
    
    p_ma = p(find(r>=param.r_eel,1));
    p_apv = p(find(r>=param.r_eel+param.t_a,1));

end