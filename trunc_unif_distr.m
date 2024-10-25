function F = trunc_unif_distr(range)

    a = range(1);
    b = range(2);
    c = range(3);
    
    x = linspace(a,b,c);
    Fcdf = (unifcdf(x,a,b)-unifcdf(a,a,b))/(unifcdf(b,a,b)-unifcdf(a,a,b));
    F = makedist('PiecewiseLinear', 'x', x, 'Fx', Fcdf);

end