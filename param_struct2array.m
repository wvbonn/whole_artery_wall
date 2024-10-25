% convert parameter structure into array suitable for boundary-value
% problem solver. convert lengths into mm, pressures into mmHg

function p = param_struct2array(v,p0)

    p(1) = 1e6*v.k_i;
    p(2) = 1e6*v.k_a;
    p(3) = v.t_p;
    p(4) = v.p_lu;
    p(5) = 10*v.l_pv;
    p(6) = v.n_v;
    p(7) = v.n_l;
    p(8) = 1e6*v.q_l;

    p(9) = p0.t_i;
    p(10) = p0.t_a;
    p(11) = p0.mu/133;
    p(12) = 1e6*p0.k_p;
    p(13) = p0.p_v-p0.sigma*(p0.pi_v-p0.pi_ref); % transcapillary pressure difference at p = 0 in tissue
    p(14) = p0.d_v*pi;
    p(15) = p0.s; % mesh size

    if isfield(v,'R_d') && isfield(p0,'D')
        p(16) = v.R_d; 
        p(17) = p0.D; 
        p(18) = p(16)*p(17)/(p(6)*p(14)*(p(10)+p(3))^2); % kappa_vv: diffusivity across vasa vasorum endothelium
    end

end