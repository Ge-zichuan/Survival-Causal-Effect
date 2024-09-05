rng('shuffle');
%%% basic information of numerical study
filenum = 167; % ID to track simulations.
simulations = 2;
samplesize = 2000;
q = 9;
d = 1;
p = 400;
iniTune = 0.1;
groupNumbering = [0 1]; % number of groups.
options = {num2str(p);num2str(q);num2str(d);num2str(samplesize);'Jiang9New';'ma91';'AFT2';'coxt';'gamm';'0.2';num2str(iniTune)};
%          predictor X ;factors   ;d         ;samplesize         ; F    ;beta  ;t1    ;t0    ;censoring;censoring rate
 
bands = zeros(4,5);
bands(1,3:end) = [1,15,15]; % larger, larger emp var, bias, est var; good vic
bands(2,2:end) = [40,35,4,4]; % larger, larger emp var, smaller est var
bands(3,2:end) = [20,20,12,12]; % smaller, larger emp variance, bias, smaller est variance
bands(4,2:end) = [5,5,5,5]; % larger, larger emp var, bias, est var;larger, smaller emp var, est var;

% Variables to store the data
xallsave = zeros(simulations,samplesize,p);
fallsave = zeros(simulations,samplesize,q);
trueBall = zeros(p,q);
zallsave = zeros(simulations,samplesize);
tallsave = zeros(simulations,samplesize);
cenallsave = zeros(simulations,samplesize);
trtallsave = zeros(simulations,samplesize);
deltaallsave = zeros(simulations,samplesize);
FEstResult = zeros(simulations,samplesize,q);
FbadEstResult = zeros(simulations,samplesize,q);
betaEstResult = zeros(simulations,d*(q-d));
betaStdEstResult = zeros(simulations,d*(q-d));
alphaResult = zeros(simulations,q);
rhoResult = zeros(simulations,1);
rhoResultTrue = zeros(simulations,1);
rhoStdResult = zeros(simulations,1);
flags = zeros(simulations,2);
censorInfo = zeros(simulations,2);
vicResult = zeros(simulations,1);

% optimization algorithm options
options1 = optimoptions('fsolve','Display','off',...
    'Algorithm','trust-region-dogleg','MaxIter', 10000,'MaxFunEvals',30000,...
    'TolFun',1e-5);
options2 = optimoptions('fminunc','Display','off');


% Generate matrix B (the loading matrix)
if (exist(strcat('trueB',num2str(filenum),'.mat'))~=0)
    load(strcat('trueB',num2str(filenum),'.mat'));
else
    dataZ = mvnrnd(zeros(p,1),(0.5).^abs(repmat(0:(p-1),p,1)-repmat(0:(p-1),p,1)'),samplesize);
    [V,~] = eigs(dataZ*dataZ',q);
    trueB = (dataZ'*V/sqrt(q));
    trueB = trueB*diag(sign(trueB(1,:)));
    save(strcat('trueB',num2str(filenum),'.mat'),'trueB');
end
% Generate propensity coefficient 
if (exist(strcat('trueA',num2str(filenum),'.txt'))~=0)
    coefficientA = dlmread(strcat('trueA',num2str(filenum),'.txt'));
else
    coefficientA = randn(q,1);
    dlmwrite(strcat('trueA',num2str(filenum),'.txt'),coefficientA,'delimiter','\t');
end

tt = 1;
while (tt(1) <= simulations)
    %% data generation
    %%% generate factor analysis data

    % Generate factor f's, then standardize it
    dataF = gendata_factor(options{5},samplesize);
    R = chol(cov(dataF));
    dataF = (dataF-repmat(mean(dataF,1),samplesize,1))/R;

    % Generate error u
    erroru = (mvnrnd(zeros(p,1),eye(p)/2/samplesize,samplesize));

    % Generate variable X from f, B, and u
    dataX = (dataF*trueB')+erroru;

    %%% Generate treatment data A
    [treatmentA,prob] = gendata_treatment(dataF,coefficientA);

    %%% Generate survival data
    betaLower                       = gendata_beta(q,d,options{6});
    [tall1,~]                       = gendata_t(dataF(treatmentA==1,:),sum(treatmentA==1),options{7},[eye(d);betaLower]);
    [cenall1,deltaall1,kval,~] 	= gendata_cen(tall1,dataF(treatmentA==1,:),options{9},str2num(options{10}),sum(treatmentA==1));
    [tall2,~]                       = gendata_t(dataF(treatmentA==0,:),sum(treatmentA==0),options{8},[eye(d);betaLower]);
    [cenall2,deltaall2,kval2,~] = gendata_cen(tall2,dataF(treatmentA==0,:),options{9},str2num(options{10}),sum(treatmentA==0));
    tall = [tall1;tall2];
    zall = tall;
    cenall = tall;
    deltaall = tall;
    cenall(treatmentA==1) = cenall1;cenall(treatmentA==0) = cenall2;
    tall(treatmentA==1) = tall1;tall(treatmentA==0) = tall2;
    deltaall(treatmentA==1) = deltaall1;deltaall(treatmentA==0) = deltaall2;
    zall(treatmentA==1) = min(tall1,cenall1);
    zall(treatmentA==0) = min(tall2,cenall2);

    %% Step: Factor Analysis
    dimOfFA = q;

    % Direct solution
    [estF,eigd] = eig(dataX*dataX');
    estF = sqrt(samplesize)*estF(:,samplesize:-1:(samplesize-(q-1))); % Unrotate the factor scores
    estFbad = estF;
    estF = estF*diag(sign(estF(1,:)).*sign(dataF(1,:)));
    estB = dataX'*estF/samplesize; 

    %% Step: Logistic Regression on Factors
    [alphaEst,~,statsLR] = glmfit(estF,treatmentA,'binomial','link','logit','constant','off');
    estE = logsig(estF*alphaEst);

    %% Step: Semiparametric Survival (from Ge, Ma, and Lu 2021).
    % Bandwidth selection
    ntemp = samplesize^(-1/4);
    bandez1 = mean(std(estF))*ntemp;
    ntemp = sum(treatmentA==1)^(-1/4);
    bandb = std(zall(treatmentA==1))*ntemp;
    bands(1,1) = bandb; % tuning the bandwidth in estimating \Z
    ntemp = sum(treatmentA==0)^(-1/4);
    bandb = std(zall(treatmentA==0))*ntemp;
    bands(1,2) = bandb; % tuning the bandwidth in estimating \Z
    bands(2,1) = bandez1; % tuning the bandwidth in estimating \lambda, denominator
    bands(3,1) = bandez1; % tuning the bandwidth in estimating \lambda, numerator
    bands(4,1) = bandez1; % tuning the bandwidth in estimating E(YX), numerator, denominator
    % Calculate all z_i-z_j
    bandadjb = bands(1,3);
    zdiffall = zall-zall';% {i,j} = Z_i-Z_j
    aallmat = treatmentA*ones(1,samplesize);
    kbz = kernel(zdiffall/bandb/bandadjb,1,'norm')/bandb/bandadjb;% {i,j} = Z_i-Z_j

    goodIndex = {find(treatmentA==1),find(treatmentA==0)};
    subgoodIndex1 = goodIndex{1};
    samplesize1 = length(subgoodIndex1);
    bandb = bands(1,1);
    bandadjb1 = bands(1,4);
    zdiffall1 = zall(subgoodIndex1)-zall(subgoodIndex1)';% {i,j} = Z_i-Z_j
    kbz1 = kernel(zdiffall1/bandb/bandadjb1,1,'Epanechnikov')/bandb/bandadjb1;% {i,j} = Z_i-Z_j
    subgoodIndex0 = goodIndex{2};
    samplesize0 = length(subgoodIndex0);
    bandb = bands(1,2);
    bandadjb0 = bands(1,5);
    zdiffall0 = zall(subgoodIndex0)-zall(subgoodIndex0)';% {i,j} = Z_i-Z_j
    kbz0 = kernel(zdiffall0/bandb/bandadjb0,1,'Epanechnikov')/bandb/bandadjb0;% {i,j} = Z_i-Z_j

    beta0 = betaLower+(rand((q-d),d)-0.5)*0;
    estF1 = estF(subgoodIndex1,:);estF0 = estF(subgoodIndex0,:);
    deltaall1 = deltaall(subgoodIndex1);deltaall0 = deltaall(subgoodIndex0);
    [betaEst,~,eflag,output] = fsolve(@(x)scoreBetaHDM5(x,estF1,zdiffall1,deltaall1,kbz1,...
        estF0,zdiffall0,deltaall0,kbz0,goodIndex,q,d,bands),beta0,options1);
    % Compute vic
    [vic,checkErrorAll] = calvicHDM(betaEst,estF,zall,deltaall,kbz1,kbz0,goodIndex,q,d,samplesize,bands);%q*(d)*log(samplesize)

    %% Step: Estimate hazard function and rho based on e and \beta
    btransx = estF*[eye(d);reshape(betaEst,[(q-d) d])]; % paras1 include upper d-by-d block.
    gvec = [estE btransx];

    [lambdaTilde1,~] = lambdaTildeEst(sort(zall(goodIndex{1})),gvec,gvec(goodIndex{1},:),zall(goodIndex{1}),...
        deltaall(goodIndex{1}),bands,1);% {i,t}: \bb\trans\X_i, time t
    [expectMeanSurvival1,~] = estEMSTilde(lambdaTilde1,zall(goodIndex{1}),size(lambdaTilde1,1));
    [lambdaTilde0,~] = lambdaTildeEst(sort(zall(goodIndex{2})),gvec,gvec(goodIndex{2},:),zall(goodIndex{2}),...
        deltaall(goodIndex{2}),bands,0);% {i,t}: \bb\trans\X_i, time t
    [expectMeanSurvival0,~] = estEMSTilde(lambdaTilde0,zall(goodIndex{2}),size(lambdaTilde0,1));
    estRho = mean(expectMeanSurvival1)-mean(expectMeanSurvival0);

    estRhoTrue = rhoTrueHDM(max(zall(intersect(goodIndex{1},find(treatmentA==1)))),btransx,options{7})...
        -rhoTrueHDM(max(zall(intersect(goodIndex{2},find(treatmentA==0)))),btransx,options{8});
    %% estimate standard deviation of rho
    [scoreMat1,checkErrorVec] = scoreBetaHDM5vec(betaEst,estF1,zdiffall1,deltaall1,kbz1,...
                        estF0,zdiffall0,deltaall0,kbz0,goodIndex,q,d,bands);
    [stdVal,jaccobMat1,~,checkErrorJaccob] = jaccobBetaNew(betaEst,estF,zall,treatmentA,deltaall,kbz1,kbz0,goodIndex,q,d,bands);
    tau1 = max(zall(goodIndex{1}));
    [Omega1,errorCheck1] = estOmega1(gvec,estF,alphaEst,zall,treatmentA,aallmat,deltaall,...
        samplesize,q,d,bands,1,tau1,scoreMat1,jaccobMat1,goodIndex,options{7});
    tau0 = max(zall(goodIndex{2}));
    [Omega0,errorCheck0] = estOmega1(gvec,estF,alphaEst,zall,1-treatmentA,1-aallmat,deltaall,...
        samplesize,q,d,bands,0,tau0,scoreMat1,jaccobMat1,goodIndex,options{8});

    %% Wrap up everything
    xallsave(tt(1),:,:) = dataX;
    fallsave(tt(1),:,:) = dataF;
    zallsave(tt(1),:) = zall;
    tallsave(tt(1),:) = tall;
    cenallsave(tt(1),:) = cenall;
    trtallsave(tt(1),:) = treatmentA;
    deltaallsave(tt(1),:) = deltaall;
    vicResult(tt(1),:) = vic;
    FEstResult(tt(1),:,:) = estF;
    FbadEstResult(tt(1),:,:) = estFbad;
    betaEstResult(tt(1),:) = betaEst(:);
    betaStdEstResult(tt(1),:) = stdVal;
    censorInfo(tt(1),:) = [kval kval2];
    rhoResult(tt(1),:) = estRho;
    rhoResultTrue(tt(1),1) = estRhoTrue;
    rhoStdResult(tt(1),:) = sqrt(median((Omega1-Omega0-estRho).^2));
    alphaResult(tt(1),:) = alphaEst;
    flags(tt(1),1) = eflag;
    tt = tt+1;
end

%% Store: data
dlmwrite(strcat('flag',num2str(filenum),'.txt'),flags,'delimiter','\t');
dlmwrite(strcat('bands',num2str(filenum),'.txt'),bands,'delimiter','\t');
dlmwrite(strcat('options',num2str(filenum),'.txt'),char(options),'delimiter','');
dlmwrite(strcat('trueBeta',num2str(filenum),'.txt'),betaLower(:),'delimiter','\t');

save(strcat('xall',num2str(filenum),'.mat'),'xallsave','-v7.3');
save(strcat('fall',num2str(filenum),'.mat'),'fallsave');
dlmwrite(strcat('zall',num2str(filenum),'.txt'),zallsave,'delimiter','\t');
dlmwrite(strcat('tall',num2str(filenum),'.txt'),tallsave,'delimiter','\t');
dlmwrite(strcat('cenall',num2str(filenum),'.txt'),cenallsave,'delimiter','\t');
dlmwrite(strcat('deltaall',num2str(filenum),'.txt'),deltaallsave,'delimiter','\t');
dlmwrite(strcat('trtall',num2str(filenum),'.txt'),trtallsave,'delimiter','\t');

save(strcat('resultF',num2str(filenum),'.mat'),'FEstResult');
save(strcat('resultFbad',num2str(filenum),'.mat'),'FbadEstResult');
dlmwrite(strcat('resultAlpha',num2str(filenum),'.txt'),alphaResult,'delimiter','\t');
dlmwrite(strcat('resultBeta',num2str(filenum),'.txt'),betaEstResult,'delimiter','\t');
dlmwrite(strcat('resultBetaStd',num2str(filenum),'.txt'),betaStdEstResult,'delimiter','\t');
dlmwrite(strcat('resultRho',num2str(filenum),'.txt'),rhoResult,'delimiter','\t');
dlmwrite(strcat('resultRhoTrue',num2str(filenum),'.txt'),rhoResultTrue,'delimiter','\t');
dlmwrite(strcat('resultRhoStdTrue',num2str(filenum),'.txt'),rhoStdResult,'delimiter','\t');
dlmwrite(strcat('resultVIC',num2str(filenum),'.txt'),vicResult,'delimiter','\t');
dlmwrite(strcat('censorInfo',num2str(filenum),'.txt'),censorInfo,'delimiter','\t');
