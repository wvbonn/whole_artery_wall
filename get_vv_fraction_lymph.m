function Q_vv_lv = get_vv_fraction_lymph(q_vv,q_lv,g)

    Q_vv = sum(q_vv.*g);
    Q_lv = sum(q_lv.*g);

    Q_vv_lv = Q_vv/Q_lv;

end