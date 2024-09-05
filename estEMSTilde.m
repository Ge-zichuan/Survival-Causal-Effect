function [expectMeanSurvival,estHazardsub]= estEMSTilde(lambda,zall,samplesize)
% estimation of expected mean survival
[sortedZ,~] = sort(zall);
timeInterval = [sortedZ(1)*ones(samplesize,1) repmat(diff(sortedZ)',samplesize,1)];
% nonZeroIndex = lambda(1,:)~=0;
estHazardsub = [zeros(length(lambda(:,1)),1) lambda];
lambdaCumTilde = [zeros(length(estHazardsub(:,1)),1) cumsum((estHazardsub(:,1:(end-1))+estHazardsub(:,2:end)).*timeInterval/2,2)];
expectMeanSurvival = sum(timeInterval.*(exp(-lambdaCumTilde(:,1:(end-1)))+exp(-lambdaCumTilde(:,2:end)))/2,2);
end