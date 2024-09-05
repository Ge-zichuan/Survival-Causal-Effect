function expectMeanSurvival = estEMS(lambda,sortedZI,timeInterval)
% estimation of expected mean survival
nonZeroIndex = lambda(1,:)~=0;
timeInterval = timeInterval(:,nonZeroIndex);
sortedlambdaTilde = lambda(:,sortedZI);
sortedlambdaTilde = [zeros(length(lambda(:,1)),1) sortedlambdaTilde(:,nonZeroIndex)];
lambdaCumTilde = [zeros(length(sortedlambdaTilde(:,1)),1) cumsum((sortedlambdaTilde(:,1:(end-1))+sortedlambdaTilde(2,end)).*timeInterval/2,2)];
expectMeanSurvival = sum(timeInterval.*(exp(-lambdaCumTilde(:,1:(end-1)))+exp(-lambdaCumTilde(:,2:end)))/2,2);
end