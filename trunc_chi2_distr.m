function F = trunc_chi2_distr(nu,range)

    a = range(1);
    b = range(2);
    c = range(3);
    
    x = linspace(a,b,c);

    Fcdf = (chi2cdf(x,nu)-chi2cdf(a,nu))/(chi2cdf(b,nu)-chi2cdf(a,nu));
    F = makedist('PiecewiseLinear', 'x', x, 'Fx', Fcdf);

end