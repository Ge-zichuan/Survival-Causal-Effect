function lambdaVec = lambdaD1TrueHDM(inputt,inputg,d,ttype)
lambdaVec = zeros(length(inputt),d);
switch ttype
    case 'cox1' % exp(1)*\bb\X\trans
    case 'cox2' % exp(1)*\bb\X\trans
    case 'sm1h' % hazard from paper SM1
    case 'texp' % truncated exp(1)*\bb\X\trans at 100
    case 'coxt' % cox model:t*exp(\bb\X\trans), study 3
        lambdaVec = exp(inputg).*repmat(inputt,1,d);
    case 'reciprocal' % 10*1/t*exp(\bb\X\trans), study 3 cont.
        templong = sum(exp(inputg),2);
        for i = 1:d
            lambdaVec(:,i) = templong.*10./(inputt+1);
        end
    case 'xia1' % Model 1 from Xia's paper
    case 'xia2' % Model 2 from Xia's paper
    case 'xia3' % Model 3 from Xia's paper
    case 'logl' % Log logistic model of T, parameter (exp(\bb\trans\x/2),2), study 2
        templong = sum(exp(inputg),2);
        lambdaVec = exp(inputg).*repmat(2*inputt./(1+inputt.^2.*templong).^2,1,d);
    case 'add1' % additive model from Chen 2007
    case 'prp1' % proportional mean residual life model from Chen 2005a
    case 'lnl2' % Log logistic model of T, parameter (exp(\bb\trans\x/2),2), study 1
end

end