function matchedData = npMatch(parasLower,estF,estE,zall,treatmentA,deltaall,samplesize,q,d,bands)
% groupIndex: 1: transplant, 0: nontransplant
matchIndex = zeros(samplesize,2);% each column match the same column
for i = 1:samplesize
    if (treatmentA(i) == 1)
        [~,minInd] = min(estF(treatmentA==0,:)*(estF(i,:)'));
        matchIndex(i,2) = minInd;
    else
        [~,minInd] = min(estF(treatmentA==1,:)*(estF(i,:)'));
        matchIndex(i,1) = minInd;
    end
end
matchedzall1 = zall;matchedzall1(treatmentA==0) = zall(matchIndex(matchIndex(:,1)~=0,1));
matchedzall0 = zall;matchedzall0(treatmentA==1) = zall(matchIndex(matchIndex(:,2)~=0,2));
matcheddeltaall1 = deltaall;matcheddeltaall1(treatmentA==0) = deltaall(matchIndex(matchIndex(:,1)~=0,1));
matcheddeltaall0 = deltaall;matcheddeltaall0(treatmentA==1) = deltaall(matchIndex(matchIndex(:,2)~=0,2));
matchedestF1 = estF;matchedestF1(treatmentA==0,:) = estF(matchIndex(matchIndex(:,1)~=0,1),:);
matchedestF0 = estF;matchedestF0(treatmentA==1,:) = estF(matchIndex(matchIndex(:,2)~=0,2),:);
% nonparametric estimation
btransx = matchedestF1*[eye(d);reshape(parasLower,[(q-d) d])]; % paras1 include upper d-by-d block.
gvec = [matchedestF1];
lambdaTilde1 = lambdaTildeEst(sort(matchedzall1),gvec,gvec,matchedzall1,...
    matcheddeltaall1,bands,1);% {i,t}: \bb\trans\X_i, time t
[expectMeanSurvival1,~] = estEMSTilde(lambdaTilde1,matchedzall1,samplesize);
btransx = matchedestF0*[eye(d);reshape(parasLower,[(q-d) d])]; % paras1 include upper d-by-d block.
gvec = [matchedestF0];
lambdaTilde0 = lambdaTildeEst(sort(matchedzall0),gvec,gvec,matchedzall0,...
    matcheddeltaall0,bands,0);% {i,t}: \bb\trans\X_i, time t
[expectMeanSurvival0,~] = estEMSTilde(lambdaTilde0,matchedzall0,samplesize);
matchedData = mean(expectMeanSurvival1)-mean(expectMeanSurvival0);

end