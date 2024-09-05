function f = funAFT1long(x,b,constant)
    f = 1-normcdf(log(x/100)-b,0,constant);
end