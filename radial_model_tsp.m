% read param_struct2array.m for correspondence between array index (p) and
% parameter names

function [c,dc] = radial_model_tsp(v,q_vv,q_lv,varpar,param) % here, p is parameter array, y is pressure

    p = param_struct2array(varpar,param);
    % calculate fluid-velocity gradient
    n = numel(varpar.r);
    dv(2:n) = (v(2:n)-v(1:n-1))/p(15);
    dv(1) = dv(2);
    assignin('base','dv',dv);
   
    sol_i = bvpinit(rmesh(p),@(r)yinit(r,p));
    sol = bvp4c(@(r,y)f(r,y,varpar.r,v,dv,q_vv,q_lv,p), @bc, sol_i); 
    
    % interpolate on regular grid
    c = interp1(sol.x,sol.y(1,:),varpar.r);
    dc = interp1(sol.x,sol.y(2,:),varpar.r);
    
end


% ADR equation with vv and lv sinks
function dydr = f(r,y,r0,v,dv,q_vv,q_lv,p) 
    
    % vv convective + vv diffusive + lv convective 
    % we ensure with min([r r0(end)]) that r is not larger than the
    % external PVAT radius - machine precision may cause this to happen and
    % cause interp1 to return NaN
    dydr = [y(2)
            -1/r*y(2) + interp1(r0,v,min([r r0(end)]))/p(17)*y(2) + interp1(r0,dv,min([r r0(end)]))/p(17)*y(1) +...
             interp1(r0,q_vv,min([r r0(end)]))/p(17)*y(1) + p(18)*p(6)*p(14)/p(17)*y(1) + interp1(r0,q_lv,min([r r0(end)]))/p(17)*y(1) ]; 
    
end

% boundary conditions
function res = bc(YL,YR)

    res = [YL(1)-1;
           YR(2)];

end

% initial mesh
function r = rmesh(p)
    
    r = 2:p(15):2+p(10)+p(3);

end

% initial guess (linear concentration)
function y0 = yinit(r,p)

    dc = 1;
    y0 = [1-dc*(r-2)/(p(10)+p(3));
          -dc/(p(10)+p(3))];

end    

