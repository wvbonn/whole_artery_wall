function [dil,d_dil10] = get_dilution_distance(v,q_vv,r,param)

    % dilution coefficient of lumen-originating fluid
    m_vv = q_vv*param.s*1e3; % mass flux from vv in µm/s
    dil(1) = 1;
    for i = 2:numel(r)
        if m_vv(i) > 0 
            dil(i) = dil(i-1)*v(i)/(v(i)+m_vv(i));
        else
            dil(i) = dil(i-1);
        end
    end
    % distance at which the fluid is diluted to less than 10%
    d_dil10 = find(dil<0.1,1);
    if numel(d_dil10)==0
        d_dil10 = r(end);
    else
        d_dil10 = r(d_dil10);
    end
    
end