function y = genData_fun_coxTang1(x,etaVec,templong)
lambda0 = 4;
beta0 = 0.5;
lambdaVec = ones(length(etaVec),1);
partProd = prod(1-x*(etaVec./lambdaVec'),2);
y = partProd.*exp((templong-lambda0*exp(beta0)).*x);
end
