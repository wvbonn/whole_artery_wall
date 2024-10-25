function F = trunc_f_distr(range,ndf,ddf,ct)

    a = range(1);
    b = range(2);
    c = range(3);
    
    x = linspace(a,b,c);

    Fcdf = (ncfcdf(x,ndf,ddf,ct)-ncfcdf(a,ndf,ddf,ct))/(ncfcdf(b,ndf,ddf,ct)-ncfcdf(a,ndf,ddf,ct));
    F = makedist('PiecewiseLinear', 'x', x, 'Fx', Fcdf);

end