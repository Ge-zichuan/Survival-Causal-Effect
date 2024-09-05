function f = funAFT1(x,b,constant)
    f = 1-normcdf(log(x)-b,0,constant);
end