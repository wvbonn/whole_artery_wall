function v = get_inner_layer_velocity(v,r,param)

    id_ma = find(r>=param.r_eel,1);
    v = sum(v(1:id_ma).*r(1:id_ma))/sum(r(1:id_ma));

end