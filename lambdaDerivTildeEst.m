function lambdaDerivVec = lambdaDerivTildeEst(inputt,inputg,zall,aall,deltaall,samplesize,goodInd,bands,glength,kerD)
% This is the estimation of lambda, rows are i, columns are t.
% khbetax: K_h(\bb\trans\x)
% denoVev: \sumi I()K_h()
ng = size(inputg);
nt = length(inputt);
numVec1 = zeros(nt,length(inputg),glength);
numVec3 = zeros(nt,length(inputg),glength);
deltaallmat = repmat(deltaall,1,nt);
yMinusT = repmat(zall,1,nt)-repmat(inputt',samplesize,1);% {i,j} = Z_i - t_j
aallmat = aall*ones(1,nt);
% read bandwidth and bandwidth adjustment
banddh = bands(2,1);
bandadjdhn = bands(2,2);
bandadjdhd = bands(2,3);
bandadjdhn1 = bands(2,4);
bandadjdhd1 = bands(2,5);
bandnh = bands(3,1);
bandadjnhn = bands(3,2);
bandadjnhd = bands(3,3);
bandadjnhn1 = bands(3,4);
bandadjnhd1 = bands(3,5);
% calculate all \bb\trans\x_i-\bb\trans\x_j (i row j column)
inputgall = permute(repmat(inputg',[1 1 samplesize]),[3 2 1])-permute(repmat(inputg,[1 1 samplesize]),[1 3 2]);% {i,j} = \bb\trans\x_j-\bb\trans\x_i
% calculate nonparametric \lambda 
bandTwo = permute(repmat([banddh*bandadjdhd banddh*bandadjdhd1]',[1 ng(1) samplesize]),[3 2 1]);
khbetax = kernel(inputgall./bandTwo,3,'Epanechnikov')/(banddh^ng(2)*bandadjdhd*bandadjdhd1);
denoVec = khbetax(:,goodInd)*(aallmat(goodInd,:).*(yMinusT(goodInd,:)>=0));
bandTwo = permute(repmat([banddh*bandadjnhd banddh*bandadjnhd1]',[1 ng(1) samplesize]),[3 2 1]);
khbetax = kernel(inputgall./bandTwo,3,'Epanechnikov')/(bandnh^ng(2)*bandadjnhd*bandadjnhd1);
khprimebetax = kernel2(inputgall/bandnh/bandadjnhn,kerD,ng(2),'Epanechnikov')/(bandnh*bandadjnhn)^(ng(2)+1);

for i = 1:glength
    numVec1(:,:,i) = (khprimebetax(:,goodInd,i)*(aallmat(goodInd,:).*(yMinusT(goodInd,:)==0).*deltaallmat(goodInd,:)))';
    numVec3(:,:,i) = (khprimebetax(:,goodInd,i)*(aallmat(goodInd,:).*(yMinusT(goodInd,:)>=0)))';
end
numVec2 = (khbetax(:,goodInd)*(aallmat(goodInd,:).*(yMinusT(goodInd,:)==0).*deltaallmat(goodInd,:)))';
lambdaVec = numVec1./repmat(denoVec',[1 1 glength])-repmat(numVec2,[1 1 glength]).*numVec3./(repmat(denoVec',[1 1 glength]).^2);

lambdaDerivVec = lambdaVec.*permute(repmat(inputg,[1 1 nt]),[3 1 2]);
lambdaDerivVec(isnan(lambdaDerivVec)) = 0;
end