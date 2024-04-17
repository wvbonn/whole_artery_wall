function g = get_solution_grid(v,p)

    g = p.r_eel-p.t_im:p.s:p.r_eel+p.t_a+v.t_pv;

end