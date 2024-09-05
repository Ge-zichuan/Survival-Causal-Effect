function lambdaVec = lambdaEstHDMsimple2(inputt,inputg,obst,obsg,aall,deltaall,samplesize,d,bands,a)
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
bandadjdhn1 = bands(2,4);
bandadjdhd1 = bands(2,5);
% calculate all Z_i-estZ
% aallmat = aall*ones(1,length(inputt));
zdiffall = obst*ones(1,length(obst))-ones(length(obst),1)*obst';% {i,j} = Z_i-Z_j
zdiffest = obst*ones(1,length(inputt))-ones(length(obst),1)*inputt';% {i,j} = Z_i-Z_j
kbz = kernel(zdiffest/(bandb*bandadjb),1,'Epanechnikov')/(bandb*bandadjb);% {i,j} = Z_i-Z_j
% calculate all \bb\trans\x_i-\bb\trans\x_j (i row j column)
if (d>1)
    kerD = 3;
    btxall = permute(repmat(inputg',[1 1 length(obst)]),[3 2 1])-permute(repmat(obsg,[1 1 estSample]),[1 3 2]);% {i,j} = \bb\trans\x_j-\bb\trans\x_i
else
    kerD = 1;
    btxall = repmat(obsg,[1 estSample])-repmat(inputg',[length(obst) 1]);% {i,j} = \bb\trans\x_j-\bb\trans\x_i
end
bandTwo = permute(repmat([banddh*bandadjdhd banddh*bandadjdhd1]',[1 estSample length(obst)]),[3 2 1]);
khbetax = kernel(btxall./bandTwo,kerD,'Epanechnikov')/(banddh^d*bandadjdhd*bandadjdhd1);
denoVec = (zdiffall<=0)*khbetax;% j row, i column
bandTwo = permute(repmat([banddh*bandadjdhn banddh*bandadjdhn1]',[1 estSample length(obst)]),[3 2 1]);
khbetax = kernel(btxall./bandTwo,kerD,'Epanechnikov')/(banddh^d*bandadjdhn*bandadjdhn1);
lambdaVec = (kbz'*(repmat(deltaall,1,estSample).*khbetax./denoVec))';
end