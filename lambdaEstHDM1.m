function lambdaVec = lambdaEstHDM1(inputt,inputg,btransx,zall,aall,deltaall,samplesize,d,bands,a)
estSample = length(inputg(:,1));
% read bandwidth and bandwidth adjustment
bandb = bands(1,1);
if a == 1
    bandadjb = bands(1,4);
elseif a == 0
    bandadjb = bands(1,5);
end
banddh = bands(2,1);
bandadjdhn = bands(2,2);
bandadjdhd = bands(2,3);
% calculate all Z_i-estZ
zdiffall = zall*ones(1,samplesize)-ones(samplesize,1)*zall';% {i,j} = Z_i-Z_j
zdiffest = zall*ones(1,estSample)-ones(samplesize,1)*inputt';% {i,j} = Z_i-Z_j
kbz = kernel(zdiffest/bandb/bandadjb,1,'norm')/bandb/bandadjb;% {i,j} = Z_i-Z_j
% calculate \bb\trans\x
% calculate all \bb\trans\x_i-\bb\trans\x_j (i row j column)
if (d>1)
    btxall = permute(repmat(inputg',[1 1 samplesize]),[3 2 1])-permute(repmat(btransx,[1 1 estSample]),[1 3 2]);% {i,j} = \bb\trans\x_j-\bb\trans\x_i
else
    btxall = repmat(btransx,[1 estSample])-repmat(inputg',[samplesize 1]);% {i,j} = \bb\trans\x_j-\bb\trans\x_i
end
% calculate estimating equation, each term.

% calculate nonparametric \lambda and \lambda\prime
if (d>1)
    kerD = 3;
else
    kerD = 1;
end
khbetax = kernel(btxall/banddh/bandadjdhd,kerD,'norm')/(banddh*bandadjdhd)^d;
denoVec = (zdiffall<=0)*khbetax;
% denoVec = denoVec(goodI,:);
khbetax = kernel(btxall/banddh/bandadjdhn,kerD,'norm')/(banddh*bandadjdhn)^d;
% lambdaVec = sum(kbz(goodI,:).*repmat(deltaall(goodI).*aall(goodI),1,estSample).*khbetax(goodI,:)./denoVec,1)';
lambdaVec = sum(kbz.*repmat(deltaall.*aall,1,estSample).*khbetax./denoVec,1)';
end