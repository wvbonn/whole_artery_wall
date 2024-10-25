function u = get_inner_layer_velocity(u,r,param)

    id_ia = find(r>=param.r_ia,1);
    u = sum(u(1:id_ia).*r(1:id_ia))/sum(r(1:id_ia));

end