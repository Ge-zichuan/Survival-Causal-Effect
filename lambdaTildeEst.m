function [lambdaVec, errorCheck] = lambdaTildeEst(inputt,inputg,obsg,obst,deltaall,bands,group)
% This is the estimation of lambda, rows are i, columns are t.
% khbetax: K_h(\bb\trans\x)
% denoVev: \sumi I()K_h()
errorCheck = strings(4,2);
ng = size(inputg);
nt = length(inputt);
deltaallmat = repmat(deltaall,1,nt);
yMinusT = obst-inputt';% {i,j} = Z_i - t_j
% aallmatinner = aall*ones(1,nt);
% read bandwidth and bandwidth adjustment
banddh = bands(2,1);
if group == 1
    bandadjdhn = bands(2,2);
    bandadjdhd = bands(2,3);
    bandadjdhn1 = bands(2,4);
    bandadjdhd1 = bands(2,5);
elseif group == 0
    bandadjdhn = bands(3,2);
    bandadjdhd = bands(3,3);
    bandadjdhn1 = bands(3,4);
    bandadjdhd1 = bands(3,5);
end
% calculate all \bb\trans\x_i-\bb\trans\x_j (i row j column)
inputgall = permute(repmat(inputg',[1 1 length(obst)]),[3 2 1])-permute(repmat(obsg,[1 1 ng(1)]),[1 3 2]);% {i,j} = \bb\trans\x_j-\bb\trans\x_i
% calculate nonparametric \lambda 
bandTwo = permute(repmat([banddh*bandadjdhd banddh*bandadjdhd1]',[1 ng(1) length(obst)]),[3 2 1]);
% inputgall./bandTwo;aall = ans;for ii = 1:size(ans,3);a = aall(:,:,ii);hist(a(:));s = waitforbuttonpress;end;
khbetax = kernel(inputgall./bandTwo,3,'norm')/(banddh^(ng(2))*bandadjdhd*bandadjdhd1);%Epanechnikov
denoVec = khbetax'*(yMinusT>=0);
if sum(isnan(denoVec),"all")>0
    errorCheck(1,1) = "denomenator (band 2,3): NaN";
    errorCheck(1,2) = 1;
else
    errorCheck(1,1) = "denomenator (band 2,3): Ok!";
    errorCheck(1,2) = 0;
end
if sum((denoVec==0),"all")>0
    errorCheck(2,1) = "denomenator (band 2,3): 0";
    errorCheck(2,2) = 1;
else
    errorCheck(2,1) = "denomenator (band 2,3): Ok!";
    errorCheck(2,2) = 0;
end
bandTwo = permute(repmat([banddh*bandadjdhn banddh*bandadjdhn1]',[1 ng(1) length(obst)]),[3 2 1]);
khbetax = kernel(inputgall./bandTwo,3,'norm')/(banddh^(ng(2))*bandadjdhn*bandadjdhn1);
if sum(isnan(khbetax),"all")>0
    errorCheck(3,1) = "numerator (band 2,2): NaN";
    errorCheck(3,2) = 1;
else
    errorCheck(3,1) = "numerator (band 2,2): Ok!";
    errorCheck(3,2) = 0;
end

numVec = khbetax'*((yMinusT==0).*deltaallmat);
lambdaVec = numVec./denoVec;
if sum(isnan(lambdaVec),"all")>0
    errorCheck(4,1) = "The ratio: NaN";
    errorCheck(4,2) = 1;
else
    errorCheck(4,1) = "The ratio: Ok!";
    errorCheck(4,2) = 0;
end
% columns are time t, rows are covariates g.

end