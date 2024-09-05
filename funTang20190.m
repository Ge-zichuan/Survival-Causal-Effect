function f = funTang20190(x,b,constant)
lambda0 = 4;
beta0 = 0;
lambdaVec = 1;
partProd = (1-x*(constant./lambdaVec)).^8;
f = partProd.*exp((-b-lambda0*exp(beta0)).*x);
end