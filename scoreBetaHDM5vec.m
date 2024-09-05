function [scoreVal,checkErrorVec] = scoreBetaHDM5vec(parasLower,xall1,zdiffall1,deltaall1,kbz1,...
    xall0,zdiffall0,deltaall0,kbz0,goodI,q,d,bands)
% This version of score function includes the treatment variable A, seperately.
% denoVev: \sumi I()K_h()
goodIndex1 = goodI{1};
goodIndex0 = goodI{2};
scoreVal = zeros(length(goodIndex1)+length(goodIndex0),(q-d)*d);
% read bandwidth and bandwidth adjustment
% calculate all \bb\trans\x_i-\bb\trans\x_j (i row j column)
% calculate nonparametric \lambda and \lambda\prime
if (d>1)
    kerD = 3;
else
    kerD = 1;
end
% subsections of score function
samplesize1 = length(goodIndex1);
btransx = xall1*[eye(d);parasLower];
if (d>1)
    btxall = permute(repmat(btransx',[1 1 samplesize1]),[3 2 1])-...
        permute(repmat(btransx,[1 1 samplesize1]),[1 3 2]);% {i,j} = \bb\trans\x_j-\bb\trans\x_i
else
    btxall = repmat(btransx,[1 samplesize1])-...
        repmat(btransx',[samplesize1 1]);% {i,j} = \bb\trans\x_j-\bb\trans\x_i
end
[scoreValallsub1,checkError1] = scoreBetaHDMsub(btxall,xall1,zdiffall1,kbz1,deltaall1,samplesize1,q,d,bands,kerD);
scoreVal(goodIndex1,:) = scoreValallsub1;

samplesize0 = length(goodIndex0);
btransx = xall0*[eye(d);parasLower];
if (d>1)
    btxall = permute(repmat(btransx',[1 1 samplesize0]),[3 2 1])-...
        permute(repmat(btransx,[1 1 samplesize0]),[1 3 2]);% {i,j} = \bb\trans\x_j-\bb\trans\x_i
else
    btxall = repmat(btransx,[1 samplesize0])-...
        repmat(btransx',[samplesize0 1]);% {i,j} = \bb\trans\x_j-\bb\trans\x_i
end
[scoreValallsub0,checkError0] = scoreBetaHDMsub(btxall,xall0,zdiffall0,kbz0,deltaall0,samplesize0,q,d,bands,kerD);
scoreVal(goodIndex0,:) = scoreValallsub0;
% scoreVal = (scoreVal+sum(scoreValallsub0,1));
checkErrorVec = [checkError1,checkError0(:,2)];
end