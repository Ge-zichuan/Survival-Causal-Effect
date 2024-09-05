function [stdVal,jaccobVal,scoreValAll,checkErrorAll] = jaccobBetaNew(parasLower,xall,zall,aall,deltaall,kbz1,kbz0,goodI,q,d,bands)
% khbetax: K_h(\bb\trans\x)
% denoVev: \sumi I()K_h()
paran = d*(q-d);
jaccobVal = zeros(length(parasLower(:)),length(parasLower(:)));
scoreValAll = zeros(length(zall),length(parasLower(:)));

goodIndex1 = (goodI{1});
goodIndex0 = (goodI{2});
% read bandwidth and bandwidth adjustment
% calculate all \bb\trans\x_i-\bb\trans\x_j (i row j column)
% calculate nonparametric \lambda and \lambda\prime
if (d>1)
    kerD = 3;
else
    kerD = 1;
end
% subsections of score function
subgoodIndex1 = intersect(goodIndex1,find(aall==1));
samplesize1 = length(subgoodIndex1);
zdiffall = zall(subgoodIndex1)*ones(1,samplesize1)-ones(samplesize1,1)*zall(subgoodIndex1)';% {i,j} = Z_i-Z_j
if (d>1)
    btxall = permute(repmat((xall(subgoodIndex1,:)*[eye(d);parasLower])',[1 1 samplesize1]),[3 2 1])-...
        permute(repmat(xall(subgoodIndex1,:)*[eye(d);parasLower],[1 1 samplesize1]),[1 3 2]);% {i,j} = \bb\trans\x_j-\bb\trans\x_i
else
    btxall = repmat(xall(subgoodIndex1,:)*[eye(d);parasLower],[1 samplesize1])-...
        repmat((xall(subgoodIndex1,:)*[eye(d);parasLower])',[samplesize1 1]);% {i,j} = \bb\trans\x_j-\bb\trans\x_i
end
[scoreValallsub1,checkError1] = scoreBetaHDMsub(btxall,xall(subgoodIndex1,:),zdiffall,kbz1,deltaall(subgoodIndex1),samplesize1,q,d,bands,kerD);
scoreValAll(subgoodIndex1,:) = scoreValallsub1;

subgoodIndex0 = intersect(goodIndex0,find(aall==0));
samplesize0 = length(subgoodIndex0);
zdiffall = zall(subgoodIndex0)*ones(1,samplesize0)-ones(samplesize0,1)*zall(subgoodIndex0)';% {i,j} = Z_i-Z_j
if (d>1)
    btxall = permute(repmat((xall(subgoodIndex0,:)*[eye(d);parasLower])',[1 1 samplesize0]),[3 2 1])-...
        permute(repmat(xall(subgoodIndex0,:)*[eye(d);parasLower],[1 1 samplesize0]),[1 3 2]);% {i,j} = \bb\trans\x_j-\bb\trans\x_i
else
    btxall = repmat(xall(subgoodIndex0,:)*[eye(d);parasLower],[1 samplesize0])-...
        repmat((xall(subgoodIndex0,:)*[eye(d);parasLower])',[samplesize0 1]);% {i,j} = \bb\trans\x_j-\bb\trans\x_i
end
[scoreValallsub0,checkError0] = scoreBetaHDMsub(btxall,xall(subgoodIndex0,:),zdiffall,kbz0,deltaall(subgoodIndex0),samplesize0,q,d,bands,kerD);
scoreValAll(subgoodIndex0,:) = scoreValallsub0;

for i = 1:paran
    for j = 1:paran
        jaccobVal(i,j) = mean(scoreValAll(:,i).*scoreValAll(:,j));
    end
end

jaccobVal = jaccobVal+0.000*eye(paran);
% sqrt(diag(inv(jaccobVal)*expectVal))/sqrt(samplesize)
stdVal = sqrt(diag(inv(jaccobVal))/length(scoreValAll(:,1)));

checkErrorAll = [checkError1 checkError0(:,2)];
end