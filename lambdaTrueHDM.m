function [lambda,lambdaCum,m,integralS,cumLambdaPrime,lambdadt,lambdaPrime,SC] = lambdaTrueHDM(inputt,inputg,ttype,kval)
% true lambda, Lambda, mrl, S from 0 to t, Lambda',dlambda/dt,lambda',sc
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
        constant = 0.1;
        templong = sum(exp(inputg),2);
        lambda = 2*(constant^2*inputt).*templong./(1+(constant*inputt).^2.*templong);
        lambdadt = 2*constant^2*templong.*(1-(constant*inputt).^2.*templong)./(1+(constant*inputt).^2.*templong).^2;
        lambdaCum = log(1+(constant*inputt).^2.*templong);
        lambdaPrime = exp(inputg).*repmat(2*constant^2*inputt./(1+(constant*inputt).^2.*templong).^2,1,length(inputg(1,:)));
        m = (1+(constant*inputt).^2.*templong)./sqrt(constant^2*templong).*(pi/2-atan((constant*inputt).*sqrt(templong)));
        integralS = atan(constant*inputt.*sqrt(templong))./sqrt(constant^2*templong);
        cumLambdaPrime = repmat((constant*inputt).^2./(1+(constant*inputt).^2.*templong),1,length(inputg(1,:))).*exp(inputg);
        SC = 1-gamcdf(inputt,15,kval);
    case 'coxtTang0' % Log logistic model of T, parameter (exp(\bb\trans\x/2),2)
        constant = -1;
        templong = inputg;
        lambda = 0;
        lambdadt = 0;
        lambdaCum = 0;
        lambdaPrime = 0;
        m = 0;
        integralS = inputt;
        for i = 1:length(inputt)
            integralS(i) = integral(@(x)funTang20190(x,templong(i),constant),0,inputt(i));
        end
        cumLambdaPrime = 0;
        SC = 1-gamcdf(inputt,15,kval);
    case 'coxtTang1' % Log logistic model of T, parameter (exp(\bb\trans\x/2),2)
        constant = -1;
        templong = inputg;
        lambda = 0;
        lambdadt = 0;
        lambdaCum = 0;
        lambdaPrime = 0;
        m = 0;
        integralS = inputt;
        for i = 1:length(inputt)
            integralS(i) = integral(@(x)funTang20191(x,templong(i),constant),0,inputt(i));
        end
        cumLambdaPrime = 0;
        SC = 1-gamcdf(inputt,15,kval);
    case 'add1' % additive model from Chen 2007
    case 'prp1' % proportional mean residual life model from Chen 2005a
    case 'lnl2' % Log logistic model of T, parameter (exp(\bb\trans\x/2),2), study 1
    case 'AFT1' % proportional mean residual life model from Chen 2005a
        constant = 0.1;
        templong = sum(logsig(inputg),2)-0.5;
        lambda = normpdf(log(inputt)-templong,0,constant)./inputt./(1-normcdf(log(inputt)-templong,0,constant));
        lambdadt = 2*constant^2*templong.*(1-(constant*inputt).^2.*templong)./(1+(constant*inputt).^2.*templong).^2;
        lambdaCum = -log(1-normcdf(log(inputt)-templong,0,constant));
        lambdaPrime = exp(inputg).*repmat(2*constant^2*inputt./(1+(constant*inputt).^2.*templong).^2,1,length(inputg(1,:)));
        m = (1+(constant*inputt).^2.*templong)./sqrt(constant^2*templong).*(pi/2-atan((constant*inputt).*sqrt(templong)));
        integralS = inputt;
        for i = 1:length(inputt)
            integralS(i) = integral(@(x)funAFT1(x,templong(i),constant),0,inputt(i));
        end
        cumLambdaPrime = repmat((constant*inputt).^2./(1+(constant*inputt).^2.*templong),1,length(inputg(1,:))).*exp(inputg);
        SC = 1-gamcdf(inputt,15,kval);
    case 'AFT1long' % proportional mean residual life model from Chen 2005a
        constant = 0.1;
        templong = sum(logsig(inputg),2)-0.5;
        lambda = normpdf(log(inputt/100)-templong,0,constant)./inputt./(1-normcdf(log(inputt/100)-templong,0,constant));
        lambdadt = 2*constant^2*templong.*(1-(constant*inputt).^2.*templong)./(1+(constant*inputt).^2.*templong).^2;
        lambdaCum = -log(1-normcdf(log(inputt/100)-templong,0,constant));
        lambdaPrime = exp(inputg).*repmat(2*constant^2*inputt./(1+(constant*inputt).^2.*templong).^2,1,length(inputg(1,:)));
        m = (1+(constant*inputt).^2.*templong)./sqrt(constant^2*templong).*(pi/2-atan((constant*inputt).*sqrt(templong)));
        integralS = inputt;
        for i = 1:length(inputt)
            integralS(i) = integral(@(x)funAFT1long(x,templong(i),constant),0,inputt(i));
        end
        cumLambdaPrime = repmat((constant*inputt).^2./(1+(constant*inputt).^2.*templong),1,length(inputg(1,:))).*exp(inputg);
        SC = 1-gamcdf(inputt,15,kval);
    case 'AFT2' % Log logistic model of T, parameter (exp(\bb\trans\x/2),2), study 1
        constant = 2;
        templong = sum(logsig(inputg),2)-0.5;
        lambda = normpdf(log(inputt)-templong,0,constant)./inputt./(1-normcdf(log(inputt)-templong,0,constant));
        lambdadt = 2*constant^2*templong.*(1-(constant*inputt).^2.*templong)./(1+(constant*inputt).^2.*templong).^2;
        lambdaCum = -log(1-normcdf(log(inputt)-templong,0,constant));
        lambdaPrime = 0;%exp(templong).*repmat(2*constant^2*inputt./(1+(constant*inputt).^2.*templong).^2,1,length(inputg(1,:)));
        m = (1+(constant*inputt).^2.*templong)./sqrt(constant^2*templong).*(pi/2-atan((constant*inputt).*sqrt(templong)));
        integralS = inputt;
        for i = 1:length(inputt)
            integralS(i) = integral(@(x)funAFT1(x,templong(i),constant),0,inputt(i));
        end
        cumLambdaPrime = repmat((constant*inputt).^2./(1+(constant*inputt).^2.*templong),1,length(inputg(1,:))).*exp(inputg);
        SC = 1-gamcdf(inputt,15,kval);
    case 'AFT4' % Log logistic model of T, parameter (exp(\bb\trans\x/2),2), study 1
        constant = 4;
        templong = sum(logsig(inputg),2)-0.5;
        lambda = normpdf(log(inputt)-templong,0,constant)./inputt./(1-normcdf(log(inputt)-templong,0,constant));
        lambdadt = 2*constant^2*templong.*(1-(constant*inputt).^2.*templong)./(1+(constant*inputt).^2.*templong).^2;
        lambdaCum = -log(1-normcdf(log(inputt)-templong,0,constant));
        lambdaPrime = exp(inputg).*repmat(2*constant^2*inputt./(1+(constant*inputt).^2.*templong).^2,1,length(inputg(1,:)));
        m = (1+(constant*inputt).^2.*templong)./sqrt(constant^2*templong).*(pi/2-atan((constant*inputt).*sqrt(templong)));
        integralS = inputt;
        for i = 1:length(inputt)
            integralS(i) = integral(@(x)funAFT1(x,templong(i),constant),0,inputt(i));
        end
        cumLambdaPrime = repmat((constant*inputt).^2./(1+(constant*inputt).^2.*templong),1,length(inputg(1,:))).*exp(inputg);
        SC = 1-gamcdf(inputt,15,kval);
    case 'AFT4long' % Log logistic model of T, parameter (exp(\bb\trans\x/2),2), study 1
        constant = 4;
        templong = sum(logsig(inputg),2)-0.5;
        lambda = normpdf(log(inputt/100)-templong,0,constant)./inputt./(1-normcdf(log(inputt/100)-templong,0,constant));
        lambdadt = 2*constant^2*templong.*(1-(constant*inputt).^2.*templong)./(1+(constant*inputt).^2.*templong).^2;
        lambdaCum = -log(1-normcdf(log(inputt/100)-templong,0,constant));
        lambdaPrime = exp(inputg).*repmat(2*constant^2*inputt./(1+(constant*inputt).^2.*templong).^2,1,length(inputg(1,:)));
        m = (1+(constant*inputt).^2.*templong)./sqrt(constant^2*templong).*(pi/2-atan((constant*inputt).*sqrt(templong)));
        integralS = inputt;
        for i = 1:length(inputt)
            integralS(i) = integral(@(x)funAFT1long(x,templong(i),constant),0,inputt(i));
        end
        cumLambdaPrime = repmat((constant*inputt).^2./(1+(constant*inputt).^2.*templong),1,length(inputg(1,:))).*exp(inputg);
        SC = 1-gamcdf(inputt,15,kval);
    case 'AFT8' % Log logistic model of T, parameter (exp(\bb\trans\x/2),2), study 1
        constant = 6;
        templong = sum(logsig(inputg),2)-0.5;
        lambda = normpdf(log(inputt)-templong,0,constant)./inputt./(1-normcdf(log(inputt)-templong,0,constant));
        lambdadt = 2*constant^2*templong.*(1-(constant*inputt).^2.*templong)./(1+(constant*inputt).^2.*templong).^2;
        lambdaCum = -log(1-normcdf(log(inputt)-templong,0,constant));
        lambdaPrime = exp(inputg).*repmat(2*constant^2*inputt./(1+(constant*inputt).^2.*templong).^2,1,length(inputg(1,:)));
        m = (1+(constant*inputt).^2.*templong)./sqrt(constant^2*templong).*(pi/2-atan((constant*inputt).*sqrt(templong)));
        integralS = inputt;
        for i = 1:length(inputt)
            integralS(i) = integral(@(x)funAFT1(x,templong(i),constant),0,inputt(i));
        end
        cumLambdaPrime = repmat((constant*inputt).^2./(1+(constant*inputt).^2.*templong),1,length(inputg(1,:))).*exp(inputg);
        SC = 1-gamcdf(inputt,15,kval);
end
end