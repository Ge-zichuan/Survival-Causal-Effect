function [lambda,lambdaCum,m,integralS,cumLambdaPrime,lambdadt,lambdaPrime,SC] = lambdaTrueHDM1(inputt,inputg,ttype,kval)
% true lambda, Lambda, mrl, S from 0 to t, L',dl/dt,l',sc
% ttype---
%     cox1: ~exp(1)*\bb\X\trans;
%     sm1h: hazard from paper SM1;
%     texp: truncated exp(1)*\bb\X\trans at 100;
%     coxt: cox model with t*exp(\bb\X\trans);
%     xia1: Model 1 from Xia's paper;
%     xia2: Model 2 from Xia's paper;
%     xia3: Model 3 from Xia's paper;
switch ttype
    case 'cox1' % exp(1)*\bb\X\trans
    case 'cox2' % exp(1)*\bb\X\trans
    case 'sm1h' % hazard from paper SM1
    case 'texp' % truncated exp(1)*\bb\X\trans at 100
    case 'coxt' % cox model:t*exp(\bb\X\trans, study 3
        templong = sum(exp(inputg/length(inputg(1,:))),2);
        lambda = templong.*inputt;
        lambdadt = templong;
        lambdaCum = templong.*inputt.^2/2;
        lambdaPrime = exp(inputg/length(inputg(1,:)))/length(inputg(1,:)).*repmat(inputt,1,length(inputg(1,:)));
        m = exp(inputt.^2/2.*templong).*normcdf(-inputt,0,sqrt(1./templong)).*...
            sqrt(2*pi).*sqrt(1./templong);
        integralS = sqrt(2*pi)./sqrt(templong).*(normcdf(inputt,0,sqrt(1./templong))-0.5);
        cumLambdaPrime = exp(inputg/length(inputg(1,:)))/length(inputg(1,:)).*repmat(inputt,1,length(inputg(1,:))).^2/2;
        SC = 1-gamcdf(inputt,15,kval);
    case 'reciprocal' % 3*1/(t+1)*exp(\bb\X\trans), study 3 cont.
        const = 10;
        lambda = (const*sum(exp(inputg),2)+1)./(inputt+1);
    case 'xia1' % Model 1 from Xia's paper
        lambda = (norminv(inputt)+10)./sum(exp(inputg),2)./normpdf(norminv(inputt));
    case 'xia2' % Model 2 from Xia's paper
    case 'xia3' % Model 3 from Xia's paper
    case 'logl' % Log logistic model of T, parameter (exp(\bb\trans\x/2),2)
        templong = sum(exp(inputg),2);
        lambda = 2*inputt.*templong./(1+inputt.^2.*templong);
        lambdadt = 2*templong.*(1-inputt.^2.*templong)./(1+inputt.^2.*templong).^2;
        lambdaCum = log(1+inputt.^2.*templong);
        lambdaPrime = exp(inputg).*repmat(2*inputt./(1+inputt.^2.*templong).^2,1,length(inputg(1,:)));
        m = (1+inputt.^2.*templong)./sqrt(templong).*(pi/2-atan(inputt.*sqrt(templong)));
        integralS = atan(inputt.*sqrt(templong))./sqrt(templong);
        cumLambdaPrime = repmat(inputt.^2./(1+inputt.^2.*templong),1,length(inputg(1,:))).*exp(inputg);
        SC = 1-gamcdf(inputt,15,kval);
    case 'loglLarget' % Log logistic model of T, parameter (exp(\bb\trans\x/2),2)
        a = 0.1;
        templong = sum(exp(inputg),2);
        lambda = 2*a*inputt.*templong./(1+a*inputt.^2.*templong);
        lambdadt = 2*a*templong.*(1-a*inputt.^2.*templong)./(1+a*inputt.^2.*templong).^2;
        lambdaCum = log(1+a*inputt.^2.*templong);
        lambdaPrime = exp(inputg).*repmat(2*a*inputt./(1+a*inputt.^2.*templong).^2,1,length(inputg(1,:)));
        m = (1+a*inputt.^2.*templong)./sqrt(a*templong).*(pi/2-atan(inputt.*sqrt(a*templong)));
        integralS = atan(inputt.*sqrt(a*templong))./sqrt(a*templong);
        cumLambdaPrime = repmat(a*inputt.^2./(1+a*inputt.^2.*templong),1,length(inputg(1,:))).*exp(inputg);
        SC = 1-gamcdf(inputt,15,kval);
    case 'add1' % additive model from Chen 2007
    case 'prp1' % proportional mean residual life model from Chen 2005a
    case 'lnl2' % Log logistic model of T, parameter (exp(\bb\trans\x/2),2), study 1
end
end