function [p_ia,p_ap] = get_interface_pressures(p,r,param)
    
    p_ia = p(find(r>=param.r_ia,1));
    p_ap = p(find(r>=param.r_ia+param.t_a,1));

end