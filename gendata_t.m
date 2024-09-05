function [tout,lambdaout] = gendata_t(xinput,n,ttype,parainput)
% generate data T, n: sample size, xdim: dimension of each observation
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
        tout = -log(rand(n,1)).*sum(exp(xinput*parainput),2);
    case 'cox2' % exp(1)*\bb\X\trans
        tout = -log(rand(n,1)).*sum(exp(xinput*parainput),2);
    case 'sm1h' % hazard from paper SM1
        r1 = 0.5;
        templong = sum(exp(xinput*parainput),2);
        tout = templong.*exp(log(((1+r1)*(rand(n,1)).^(-templong*r1)-1)/r1)./templong);
    case 'texp' % truncated exp(1)*\bb\X\trans at 100
        templong = sum(exp(xinput*parainput),2);
        c = 1/(1-exp(-100*templong));
        tout = (-log(1-rand(n,1)./c)./templong);
        while (sum(tout>100) > 0)
            tout(tout>100) = -log(1-rand(sum(tout>100),1)./c)./templong(tout>100);
        end
    case 'coxt' % cox model:t*exp(\bb\X\trans, study 3
        templong = sum(exp(xinput*parainput/length(xinput(1,:))),2);
        tout = sqrt(-2*log(rand(n,1))./templong);
        lambdaout = templong.*tout;
    case 'coxtTang0' % cox model: from Tang 2019, T0
        tout = sqrt(-log(rand(n,1))/4);
        lambdaout = 4*ones(n,1);
    case 'coxtTang1' % cox model: from Tang 2019, T0
        templong = -xinput*parainput;
        etaVec = -parainput';
        options = optimoptions('fsolve','Display','off');
        tout = fsolve(@(x)genData_fun_coxTang1(x,etaVec,templong),ones(n,1),options);
        lambdaout = 4*exp(0.5)*ones(n,1);
    case 'reciprocal' % 3*1/(t+1)*exp(\bb\X\trans), study 3 cont.
        const = 10;
        templong = sum(exp(xinput*parainput),2);
        tout = rand(n,1).^(-1./(const*templong+1))-1;
        lambdaout = (const*templong+1)./(tout+1);
    case 'xia1' % Model 1 from Xia's paper
        templong = sum(exp(6*xinput*parainput+1),2);
        tout = normcdf(5*(-log(rand(n,1)).*(templong)-2));
        lambdaout = (norminv(tout)+10)./templong./normpdf(norminv(tout));
    case 'xia2' % Model 2 from Xia's paper
        templong = sum((xinput*parainput/2),2);
        tout = exp(-log(rand(n,1))+templong);
    case 'xia3' % Model 3 from Xia's paper
        templong = sum((1-xinput*parainput).^2,2);
        tout = exp(5-10*templong+normrnd(0,1,n,1));
    case 'logl' % Log logistic model of T, parameter (exp(\bb\trans\x/2),2)
        templong = sum(exp(xinput*parainput),2);
        tout = ((1./rand(n,1)-1.0)./templong).^(1/2);
%         while (sum(tout>20) > 0)
%             tout(tout>20) = templong(tout>20).*(1./rand(sum(tout>20),1)-1).^(1/8);
%         end
        lambdaout = 2*tout.*templong./(1+tout.^2.*templong);
    case 'loglLarget' % Log logistic model of T, parameter (exp(\bb\trans\x/2),2)
        templong = sum(exp(xinput*parainput),2);
        tout = ((1./rand(n,1)-1)./templong).^(1/2)/0.1;
        lambdaout = 0.02*tout.*templong./(1+(0.1*tout).^2.*templong);
    case 'add1' % additive model from Chen 2007
        templong = sum((xinput*parainput),2);
        tout = (1+templong)./sqrt(rand(n,1))-1-templong;
    case 'prp1' % proportional mean residual life model from Chen 2005a
        templong = exp(sum(xinput,2));
        tout = (rand(n,1)).^(-templong./(1+templong))-1;
    case 'lnl2' % Log logistic model of T, parameter (exp(\bb\trans\x/2),2), study 1
        templong = sum(exp(xinput*parainput),2);
        tout = ((1./rand(n,1)-1.0)./templong).^(1/2);
        lambdaout = 2*tout.*templong./(1+tout.^2.*templong);
%      case 'lnl2' % Log logistic model of T, parameter (exp(\bb\trans\x),2), study 1
%         templong = sum(exp(xinput*parainput)/2,2);
%         tout = (templong).*sqrt((1./rand(n,1)-1));
%         lambdaout = 2*tout./(templong.^2+tout.^2);
    case 'AFT1' % Log logistic model of T, parameter (exp(\bb\trans\x/2),2)
        templong = sum(logsig(xinput*parainput),2)-0.5;
        normVar = 0.1;% variation (std) of error in AFT model (normal)
        tout = exp(templong+normrnd(0,normVar,length(templong),1));
        lambdaout = normpdf(log(tout)-templong,0,normVar)./tout./(1-normcdf(log(tout)-templong,0,normVar));
    case 'AFT1long' % Log logistic model of T, parameter (exp(\bb\trans\x/2),2)
        templong = sum(logsig(xinput*parainput),2)-0.5;
        normVar = 0.1;% variation (std) of error in AFT model (normal)
        tout = exp(templong+normrnd(0,normVar,length(templong),1))*100;
        lambdaout = normpdf(log(tout/100)-templong,0,normVar)./tout./(1-normcdf(log(tout/100)-templong,0,normVar));
    case 'AFT2' % Log logistic model of T, parameter (exp(\bb\trans\x/2),2)
        templong = sum(logsig(xinput*parainput),2)-0.5;
        normVar = 2;% variation (std) of error in AFT model (normal)
        tout = exp(templong+normrnd(0,normVar,length(templong),1));
        lambdaout = normpdf(log(tout)-templong,0,normVar)./tout./(1-normcdf(log(tout)-templong,0,normVar));
    case 'AFT4' % Log logistic model of T, parameter (exp(\bb\trans\x/2),2)
        templong = sum(logsig(xinput*parainput),2)-0.5;
        normVar = 4;% variation (std) of error in AFT model (normal)
        tout = exp(templong+normrnd(0,normVar,length(templong),1));
        lambdaout = normpdf(log(tout)-templong,0,normVar)./tout./(1-normcdf(log(tout)-templong,0,normVar));
    case 'AFT4long' % Log logistic model of T, parameter (exp(\bb\trans\x/2),2)
        templong = sum(logsig(xinput*parainput),2)-0.5;
        normVar = 4;% variation (std) of error in AFT model (normal)
        tout = exp(templong+normrnd(0,normVar,length(templong),1))*100;
        lambdaout = normpdf(log(tout/100)-templong,0,normVar)./tout./(1-normcdf(log(tout/100)-templong,0,normVar));
    case 'AFT8' % Log logistic model of T, parameter (exp(\bb\trans\x/2),2)
        templong = sum(logsig(xinput*parainput),2)-0.5;
        normVar = 6;% variation (std) of error in AFT model (normal)
        tout = exp(templong+normrnd(0,normVar,length(templong),1));
        lambdaout = normpdf(log(tout)-templong,0,normVar)./tout./(1-normcdf(log(tout)-templong,0,normVar));
end
end