function [vic,checkErrorAll] = calvicHDM(paras0,xall,zall,deltaall,kbz1,kbz0,goodI,p,d,samplesize,bands)
% calculate vic
% Refer to: Yanyuan Ma, Xinyu Zhang,
% A validated information criterion to determine the structural dimension in dimension reduction models
vec1 = ones(p-1-d,1)*0.1;
vec2 = zeros(p-1-d,1);
paraExp1 = zeros(p-d-1,d+1);
paraExp2 = zeros(p-d-1,d+1);
for i = 1:(p-d-1)
    for j = 1:d
        paraExp1(i,j) = paras0(i+1,j)-vec1(i)*paras0(1,j);
        paraExp2(i,j) = paras0(i+1,j)-vec2(i)*paras0(1,j);
    end
    paraExp1(i,d+1) = vec1(i);
    paraExp2(i,d+1) = vec2(i);
end

% scoreBetaHDM4(paras0,xall,zall,deltaall,kbz1,kbz0,goodI,p,d,bands);
% scoreValnew1 = scoreBetaHDM4(paraExp1,xall,zall,deltaall,kbz1,kbz0,goodI,p-1,d+1,bands);
% scoreValnew2 = scoreBetaHDM4(paraExp2,xall,zall,deltaall,kbz1,kbz0,goodI,p-1,d+1,bands);
subgoodIndex1 = goodI{1};
subgoodIndex0 = goodI{2};
zdiffall1 = zall(subgoodIndex1)-zall(subgoodIndex1)';% {i,j} = Z_i-Z_j
zdiffall0 = zall(subgoodIndex0)-zall(subgoodIndex0)';% {i,j} = Z_i-Z_j
xall1 = xall(subgoodIndex1,:);xall0 = xall(subgoodIndex0,:);
deltaall1 = deltaall(subgoodIndex1);deltaall0 = deltaall(subgoodIndex0);
[scoreValnew1,checkerror1] = scoreBetaHDM5(paraExp1,xall1,zdiffall1,deltaall1,kbz1,...
    xall0,zdiffall0,deltaall0,kbz0,goodI,p,d+1,bands);
[scoreValnew2,checkerror2] = scoreBetaHDM5(paraExp2,xall1,zdiffall1,deltaall1,kbz1,...
    xall0,zdiffall0,deltaall0,kbz0,goodI,p,d+1,bands);
vic = (sum(scoreValnew1.^2)+sum(scoreValnew2.^2))...
    /(samplesize)/2+p*d*log(samplesize);
if isnan(vic)
    if (sum(isnan(scoreValnew1),'all')==0)
    else
        disp('scoreValnew1 has NaN');
    end
    if (sum(isnan(scoreValnew2),'all')==0)
    else
        disp('scoreValnew2 has NaN');
    end

else
end
checkErrorAll = [checkerror1 checkerror2(:,2:end)];
end

