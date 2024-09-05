function lambdaVec = cumLambdaEstMR(estt,estg,inputt,inputg,deltaall,samplesize,bands,d,group)
estgSample = length(estg(:,1));
esttSample = length(estt(:,1));
% read bandwidth and bandwidth adjustment
% if group == 1
%     bandb = bands(1,2);
%     bandadjb = bands(1,4);
% else
%     bandb = bands(1,1);
%     bandadjb = bands(1,3);
% end
banddh = bands(2,1);
bandadjdhn = bands(2,2);
bandadjdhd = bands(2,3);
% calculate all Z_i-estZ
% zdiffall = inputt*ones(1,samplesize)-ones(estSample,1)*zall';% {i,j} = Z_i-Z_j
% kbz = kernel(zdiffall/bandb/bandadjb,1,'norm')/bandb/bandadjb;% {i,j} = Z_i-Z_j
zdiffall = inputt*ones(1,samplesize)-ones(samplesize,1)*inputt';% {i,j} = Z_i-Z_j
zdiffest = inputt*ones(1,esttSample)-ones(samplesize,1)*estt';% {i,j} = Z_i-Z_j
% kbz = kernel(zdiffest/bandb/bandadjb,1,'norm')/bandb/bandadjb;% {i,j} = Z_i-Z_j
% calculate \bb\trans\x
% inputg = xall*[eye(d);parasLower]; % paras1 include upper d-by-d block.
% calculate all \bb\trans\x_i-\bb\trans\x_j (i row j column)
if (d>1)
    kerD = 3;
    btxall = permute(repmat(estg',[1 1 samplesize]),[3 2 1])-permute(repmat(inputg,[1 1 estgSample]),[1 3 2]);% {i,j} = \bb\trans\x_j-\bb\trans\x_i
else
    kerD = 1;
    btxall = repmat(inputg,[1 estgSample])-repmat(estg',[samplesize 1]);% {i,j} = \bb\trans\x_j-\bb\trans\x_i
end

% calculate nonparametric \lambda and \lambda\prime
khbetax = kernel(btxall/banddh/bandadjdhd,kerD,'norm')/banddh/bandadjdhd;
denoVec = (zdiffall<=0)*khbetax;
khbetax = kernel(btxall/banddh/bandadjdhn,kerD,'norm')/banddh/bandadjdhn;
lambdaVec = ((((zdiffest<=0)').*(repmat(deltaall,1,esttSample)'))*(khbetax./denoVec))';

end