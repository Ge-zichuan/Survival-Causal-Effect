function varRho = estVarRho(tau,estRho,xall,zdiffall,aall,aallmat,kbz,deltaall,samplesize,goodI,q,d,bands)
% khbetax: K_h(\bb\trans\x)
% denoVev: \sumi I()K_h()
evaluationT = 0:0.01:tau;
lambdaTilde = lambdaTildeEst(evaluationT,gvec,zall,aall,deltaall,samplesize,goodI{1},bands);
cumulLambdaTilde = [zeros(1,length(gvec)); cumsum((lambdaTilde(1:(end-1))+lambdaTilde(2:end))/2.*diff(evaluationT),1)];
intS1 = sum((exp(-cumulLambdaTilde(1:(end-1)))+exp(-cumulLambdaTilde(2:end)))/2.*diff(evaluationT),1);
lambdaTilde = lambdaTildeEst(evaluationT,gvec,zall,1-aall,deltaall,samplesize,goodI{2},bands);
cumulLambdaTilde = [zeros(1,length(gvec)); cumsum((lambdaTilde(1:(end-1))+lambdaTilde(2:end))/2.*diff(evaluationT),1)];
intS0 = sum((exp(-cumulLambdaTilde(1:(end-1)))+exp(-cumulLambdaTilde(2:end)))/2.*diff(evaluationT),1);

varRho = (intS1-intS0-estRho)^2;


% Omega1 = estOmega(gvec,estF,alphaEst,estbeta,xall,zall,zdiffall,aall,aallmat,kbz,deltaall,samplesize,goodI,p,q,d,bands,1);
% Omega0 = estOmega(gvec,estF,alphaEst,estbeta,xall,zall,zdiffall,aall,aallmat,kbz,deltaall,samplesize,goodI,p,q,d,bands,0);
% varRho = (Omega1-Omega0-estRho)^2;


end