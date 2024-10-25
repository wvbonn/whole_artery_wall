function F = trunc_lognorm_distr(mu,sigma,range)

    a = range(1);
    b = range(2);
    c = range(3);
    
    x = linspace(a,b,c);
    Fcdf = (logncdf(x,mu,sigma)-logncdf(a,mu,sigma))/(logncdf(b,mu,sigma)-logncdf(a,mu,sigma));
    F = makedist('PiecewiseLinear', 'x', x, 'Fx', Fcdf);

end