function Q_v_l = get_vv_fraction_lymph(q_v,q_l,g)

    Q_v = sum(q_v.*g);
    Q_l = sum(q_l.*g);

    Q_v_l = Q_v/Q_l;

end