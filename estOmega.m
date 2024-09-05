function [Omega] = estOmega(gvec,estF,alphaEst,zall,aall,aallmat...
    ,deltaall,samplesize,q,d,bands,a,tau,scoreVal,scoreSquare,goodI)
% gvec,estF,alphaEst,betaEst,zall,zdiffall,treatmentA,aallmat,kbz,deltaall,samplesize,goodIndex,p,q,d,bands,0,tau
glength = d+1;
atransf = estF*alphaEst;
% read bandwidth and bandwidth adjustment
% banddh = bands(2,1);
% bandadjdhn = bands(2,2);
% bandadjdhd = bands(2,3);
% bandnh = bands(3,1);
% bandadjnhn = bands(3,2);
% bandadjnhd = bands(3,3);
bandhe = bands(4,1);
% bandadjen = bands(4,2);
if a == 1
    bandadjed = bands(4,3);
else
    bandadjed = bands(4,5);% used for erg
end
%% Estimate: integral of e^{-cumulative hazard} from 0 to \tau
% finerCoef = 1000;
if a == 1
    modelname = 'loglLarget';
    subgoodIndex = intersect(goodI{1},find(aall==1));
else
    modelname = 'coxt';
    subgoodIndex = intersect(goodI{2},find(aall==1));
end
subsamplesize = length(subgoodIndex);
r = sum(aall)/samplesize;
[estt,esttI] = sort(zall(subgoodIndex));
[~,backI] = sort(esttI);% return to the original data
lambdaTilde = lambdaEstHDMsimple(estt,gvec,zall(subgoodIndex),gvec(subgoodIndex,:),...
    aall(subgoodIndex),deltaall(subgoodIndex),length(subgoodIndex),d+1,bands,a);
%     lambdaTilde = lambdaTildeEst(estt,gvec,gvec(subgoodIndex,:),zall(subgoodIndex),...
%         aall(subgoodIndex),deltaall(subgoodIndex),bands);% {i,t}: \bb\trans\X_i, time t
% evaluationTmat = repmat(diff(evaluationT),[1 samplesize]);
% cumulLambdaTilde = [zeros(1,samplesize); cumsum((lambdaTilde(1:(end-1),:)+lambdaTilde(2:end,:))/2.*evaluationTmat,1)];
% intS = sum((exp(-cumulLambdaTilde(1:(end-1),:))+exp(-cumulLambdaTilde(2:end,:)))/2.*evaluationTmat,1)';
% nonzeroTimeInd = lambdaTilde(1,:)~=0;
evaluationTmat = repmat(diff(estt),[1 samplesize])';
lambdaTildeTemp = lambdaTilde(:,:);
cumulLambdaTilde = [zeros(samplesize,1) cumsum((lambdaTildeTemp(:,1:(end-1))+lambdaTildeTemp(:,2:end))/2.*evaluationTmat,2)];
% intS = sum((exp(-cumulLambdaTilde(:,1:(end-1)))+exp(-cumulLambdaTilde(:,2:end)))/2.*evaluationTmat,2);
intSAll = cumLambdaEstMR(estt,gvec,zall(subgoodIndex),gvec(subgoodIndex,:),deltaall(subgoodIndex),length(subgoodIndex),bands,d+1,a);
% [~,~,~,intS,~,~,~,~] = lambdaTrueHDM(ones(samplesize,1)*tau,gvec,modelname,0.1);
% intS = integral(@(x)lambdaDouIntEstMR(parasLower,x,inputg(inputt==uniqueT(i),:),xall,zall,deltaall,samplesize,bands,d,group)...
%         ,uniqueT(i),upperLimit,'ArrayValued',true);
gveclong = repmat(gvec,[subsamplesize,1]);
timeLong = zeros(subsamplesize*samplesize,1);
for i = 1:subsamplesize
    timeLong(((i-1)*samplesize+1):(i*samplesize)) = estt(i);
end

% cumulLambdaTrue = reshape(clTrue,samplesize,length(subgoodIndex));
% lambdaTrue = reshape(lTrue,samplesize,length(subgoodIndex));
% plot(estt,lambdaTildeTemp(1,:),estt,lambdaTrue(1,:))
% [~,~,~,intStrue,~,~,~,~] = lambdaTrueHDM(ones(samplesize,1)*tau,gvec,modelname,0.1);
%% Estimate: integral of martingale multiply integral of e^{-cumulative hazard} from s to \tau
kerD = 3;
btxall = permute(repmat(gvec(subgoodIndex)',[1 1 subsamplesize]),[3 2 1])-permute(repmat(gvec(subgoodIndex),[1 1 subsamplesize]),[1 3 2]);% {i,j} = \bb\trans\x_j-\bb\trans\x_i

yMinusT = repmat(zall(subgoodIndex),1,subsamplesize)-repmat(estt',subsamplesize,1);% {i,j} = Z_i - t_j
khbetax = kernel(btxall/(bandhe*bandadjed),3,'norm')/(bandhe*bandadjed)^glength;
erg = ((yMinusT>=0)'*(khbetax))./repmat(sum(khbetax,1),[subsamplesize,1]);
erg = erg';
% reverseTemp = exp(-cumulLambdaTilde(:,end:-1:1));
% evaluationTmat = -repmat(diff(estt(end:-1:1)),[1 samplesize])';
% intSReverse = cumsum((reverseTemp(:,1:(end-1))+reverseTemp(:,2:end))/2.*evaluationTmat,2);
% intSUp = [intSReverse(:,end:-1:1) zeros(samplesize,1)];

intSUp = repmat(intS,1,subsamplesize)-[zeros(samplesize,1),cumsum((exp(-cumulLambdaTilde(:,1:(end-1)))+exp(-cumulLambdaTilde(:,2:end)))/2.*evaluationTmat,2)];
% [lTrue,clTrue,~,integralS,~,~,~,~] = lambdaTrueHDM(timeLong,gveclong(:,2),modelname,0.1);
% [~,~,~,integralSup,~,~,~,~] = lambdaTrueHDM(tau*ones(length(subgoodIndex)*samplesize,1),gveclong(:,2),modelname,0.1);
% intSUp = reshape(integralSup-integralS,samplesize,length(subgoodIndex));

indtemp = sub2ind(size(intSUp),subgoodIndex,backI);
indtemp2 = sub2ind(size(erg),(1:subsamplesize)',backI);
intMartigale1 = deltaall(subgoodIndex).*intSUp(indtemp)./erg(indtemp2);
intMartigale2inner = intSUp(subgoodIndex,:).*lambdaTilde(subgoodIndex,:)./erg;
evaluationTmat = repmat(diff(estt),[1 subsamplesize])';
intMartigale2whole = [zeros(subsamplesize,1) ...
    cumsum((intMartigale2inner(:,1:(end-1))+intMartigale2inner(:,2:end))/2.*evaluationTmat,2)];
indtemp = sub2ind(size(intMartigale2whole),(1:subsamplesize)',backI);
intMartigale2 = intMartigale2whole(indtemp);

%% Estimate: integral of e^{-cumulative hazard} multiply cumulative hazard and g_beta from 0 to \tau
if a == 1
    lambdaDerivTilde = lambdaDerivTildeEst(estt,gvec,zall,aall,deltaall,samplesize,goodI{1},bands,glength,kerD);
else
    lambdaDerivTilde = lambdaDerivTildeEst(estt,gvec,zall,aall,deltaall,samplesize,goodI{2},bands,glength,kerD);
end
elld = repmat(exp(-cumulLambdaTilde'),[1 1 glength]).*lambdaDerivTilde;
elldIntegral = squeeze(sum((elld(1:(end-1),:,:)+elld(2:end,:,:))/2.*repmat(diff(estt),[1 length(gvec(:,1)) glength]),1));
gb = cat(3,zeros(samplesize,(q-d)*d,1),repmat(repmat(estF(:,(d+1):end),[1 d]),[1,1,d]));%n x q-d x d+1
eiecgb = squeeze(sum(permute(repmat(elldIntegral,[1 1 (q-d)*d]),[1 3 2]).*gb,3));


%% Estimate: expectation of integral of e^{-cumulative hazard} multiply cumulative hazard and g_alpha from 0 to \tau
expitaf = logsig(atransf);
ga = cat(3,repmat(repmat(expitaf.*(1-expitaf),[1,q]).*estF,[1,1,1]),zeros(samplesize,q,d));
eiecga = squeeze(sum(permute(repmat(elldIntegral,[1 1 q]),[1 3 2]).*ga,3));

%% Estimate: phi
phiDenominator1 = inv(squeeze(mean(repmat(expitaf.*(1-expitaf),[1 q q]).*(repmat(estF,[1,1,q]).*permute(repmat(estF,[1,1,q]),[1 3 2])),1)));
phiDenominator = permute(repmat(phiDenominator1,[1 1 samplesize]),[3 1 2]);
phi = repmat(aall-expitaf,[1 q]).*squeeze(sum(repmat(estF,[1 1 q]).*phiDenominator,2));


%% Summarize everything
part2 = zeros(samplesize,1);
part3 = zeros(samplesize,1);
part2(subgoodIndex) = intMartigale1;
part3(subgoodIndex) = intMartigale2;
inv(scoreSquare);
Omega = intS-r^(-1)*aall.*part2-r^(-1)*aall.*part3...
    +sum(eiecgb.*(scoreVal/scoreSquare),2)...
    +sum(eiecga.*phi,2);
1;
end