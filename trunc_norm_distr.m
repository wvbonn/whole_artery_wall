function F = trunc_norm_distr(mu,sigma,range)

    a = range(1);
    b = range(2);
    c = range(3);
    
    x = linspace(a,b,c);
    Fcdf = (normcdf(x,mu,sigma)-normcdf(a,mu,sigma))/(normcdf(b,mu,sigma)-normcdf(a,mu,sigma));
    F = makedist('PiecewiseLinear', 'x', x, 'Fx', Fcdf);

end