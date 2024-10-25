function g = get_solution_grid(v,p)

    g = p.r_ia-p.t_i:p.s:p.r_ia+p.t_a+v.t_p;

end