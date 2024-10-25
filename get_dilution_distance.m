function [dil,d_dil10] = get_dilution_distance(u,q_v,r,param)

    % dilution coefficient of lumen-originating fluid
    m_v = q_v*param.s*1e3; % mass flux from vv in µm/s
    dil(1) = 1;
    for i = 2:numel(r)
        if m_v(i) > 0 
            dil(i) = dil(i-1)*u(i)/(u(i)+m_v(i));
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