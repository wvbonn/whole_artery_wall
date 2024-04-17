% read param_struct2array.m for correspondence between array index (p) and
% parameter names

function [y,v,q_vv,q_lv] = radial_model_fluid(varpar,param,ATHERO) 
% here, p is parameter array, y is pressure

    p = param_struct2array(varpar,param);
    sol_i = bvpinit(rmesh(p),@(r)yinit(r,p));
    sol = bvp4c(@(r,y)f(r,y,p,ATHERO), @(YL,YR)bc(YL,YR,p), sol_i); 
    
    % interpolate on regular grid
    y = interp1(sol.x,sol.y(1,:),varpar.r);
    v = get_velocity(sol,varpar,p);
    if ATHERO
        q_vv = p(5)*p(6)*p(14)*(p(13)-y);
    else
        q_vv = p(5)*p(6)*p(14)*(p(13)-y).*(y-p(13)<0);
    end
    q_lv = p(7)*p(8)*ones(size(y));

end


% mass conservation with vv and lv source terms
function dydr = f(r,y,p,ATHERO) 

    s = 200; % steepness of sigmoids that make permeability and microvascular fluxes C^1 functions across the whole wall
    vv_coef = p(11)*p(5)*p(6)*p(14); % mu*l_vv*A_vv
    lv_coef = p(11)*p(7)*p(8); % mu*n_lv*q_lv
    r_apv = 2+p(10);
    % continuous function for k and its 1st-order derivative
    f_k1 = p(2) + (p(1)-p(2))/(1+(r/2)^s);
    f_k = p(12) + (f_k1-p(12))/(1+(r/r_apv)^s);
    dk1_dr = -(p(1)-p(2))*s*2^-s*r^(s-1)/(1+(r/2)^s)^2;
    dk_dr = (dk1_dr*(1+(r/r_apv)^s) - ((f_k1-p(12))*s*r_apv^-s*r^(s-1))) / (1+(r/r_apv)^s)^2;

    % continuous functions for microvascular exchanges
    f_vv = vv_coef * (1-1/(1+(r/2)^s));
    f_lv = lv_coef * (1-1/(1+(r/2)^s));
    
    if ATHERO
        dydr = [y(2)
                -1/r*y(2) - dk_dr/f_k*y(2) + f_vv/f_k*(y(1)-p(13)) + f_lv/f_k ]; % disrupted glycocalyx -> Starling
    else
        dydr = [y(2)
                -1/r*y(2) - dk_dr/f_k*y(2) + f_vv/f_k*(y(1)-p(13))*(y(1)-p(13)<0) + f_lv/f_k ]; % intact glycocalyx
    end

end

% BC: Dirichlet at lumen, Neumann at external PVAT
function res = bc(YL,YR,p)

    res = [YL(1)-p(4);
           YR(2)];

end

% initial mesh
function r = rmesh(p)
    
    r = 2-p(9):p(15):2+p(10)+p(3);

end

% initial guess (linear pressure)
function y0 = yinit(r,p)

    v_ref = 3e-5; % reference fluid velocity in inner layers
    dp = p(11)*v_ref*p(9)/p(1);
    y0 = [p(4)-dp*(r-(2-p(9)))/(p(9)+p(10)+p(3));
          -dp/(p(9)+p(10)+p(3))];

end    

% get_velocity (apply Darcy's law) - output in µm/s
function v = get_velocity(sol,varpar,p)
    
    s = 200;
    g = varpar.r;
    r_apv = 2+p(10);
    f_k1 = p(2) + (p(1)-p(2))./(1+(g/2).^s);
    f_k = p(12) + (f_k1-p(12))./(1+(g/r_apv).^s);
    v = -1e3*f_k/p(11).*interp1(sol.x,sol.y(2,:),g);

end