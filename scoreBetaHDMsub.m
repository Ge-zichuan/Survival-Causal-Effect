function [scoreValall,checkError] = scoreBetaHDMsub(btxall,xall,zdiffall,kbz,deltaall,samplesize,q,d,bands,kerD)
% This version of score function includes the treatment variable A.
% khbetax: K_h(\bb\trans\x)
% denoVev: \sumi I()K_h()
checkError = strings(8,2);
scoreValall = zeros(samplesize,(q-d)*d);
% read bandwidth and bandwidth adjustment
banddh = bands(2,1);
bandadjdhn = bands(2,2);
bandadjdhd = bands(2,3);
bandnh = bands(3,1);
bandadjnhn = bands(3,2);
bandadjnhd = bands(3,3);
bandhe = bands(4,1);
bandadjen = bands(4,2);
bandadjed = bands(4,3);
% calculate nonparametric \lambda and \lambda\prime
khbetax = kernel(btxall/(banddh*bandadjdhd),kerD,'norm')/(banddh*bandadjdhd)^d;%Epanechnikov
denoVec = khbetax*(zdiffall>=0);% row i column j
if (sum(denoVec==0,'all')==0)
    checkError(1,1) = 'Deno lambda has zero (2,3):';
    checkError(1,2) = 0;
else
    checkError(1,1) = 'Deno lambda has zero (2,3):';
    checkError(1,2) = 1;
end
khbetax = kernel(btxall/(banddh*bandadjdhn),kerD,'Epanechnikov')/(banddh*bandadjdhn)^d;
% aall = ans;for ii = 1:size(ans,3);a = aall(:,:,ii);hist(a(:));s = waitforbuttonpress;end;
lambdaDen = sum(kbz.*(repmat(deltaall,1,samplesize)').*khbetax./denoVec,2);
if (sum(isnan(lambdaDen),'all')==0)
    checkError(2,1) = 'lambda has NaN:';
    checkError(2,2) = 0;
else
    checkError(2,1) = 'lambda has NaN:';
    checkError(2,2) = 1;
end
if (sum((lambdaDen==0),'all')==0)
    checkError(3,1) = 'lambda has 0 (2,2):';
    checkError(3,2) = 0;
else
    checkError(3,1) = 'lambda has 0 (2,2):';
    checkError(3,2) = 1;
end

% calculate kernel values of K\prime(.)
lambdaVec = zeros(samplesize,d);
khprimevec = kernel2(btxall/(bandnh*bandadjnhn),kerD,d,'Epanechnikov')/(bandnh*bandadjnhn)^(d+1);
khbetax = kernel(btxall/(bandnh*bandadjnhd),kerD,'Epanechnikov')/(bandnh*bandadjnhd)^d;
for t = 1:d
    lambdaVec(:,t) = (sum(-kbz.*(repmat(deltaall,1,samplesize)').*squeeze(khprimevec(:,:,t))./denoVec,2)...
        +sum(kbz.*(repmat(deltaall,1,samplesize)')./denoVec.*khbetax.*(squeeze(khprimevec(:,:,t))*(zdiffall>=0))./denoVec,2))./...
        lambdaDen;
end
if (sum(isnan(lambdaVec),'all')==0)
    checkError(4,1) = 'ratio of lambda has NaN:';
    checkError(4,2) = 0;
else
    checkError(4,1) = 'ratio of lambda has NaN:';
    checkError(4,2) = 1;
end

% calculate nonparametric expectations
khbetax = kernel(btxall/(bandhe*bandadjed),3,'Epanechnikov')/(bandhe*bandadjed)^d;
expectDen = sum(khbetax.*(zdiffall<=0),2);%./sum(khbetax,2);
if (sum(expectDen==0,'all')==0)
    checkError(5,1) = 'Deno of Expectation has 0 (4,3):';
    checkError(5,2) = 0;
else
    checkError(5,1) = 'Deno of Expectation has 0 (4,3):';
    checkError(5,2) = 1;
end

khbetax = kernel(btxall/(bandhe*bandadjen),3,'Epanechnikov')/(bandhe*bandadjen)^d;
for i = 1:d
    for j = 1:(q-d)
        scoreValall(:,(i-1)*(q-d)+j) = deltaall.*lambdaVec(:,i).*...
            (xall(:,d+j)-sum(khbetax.*repmat(xall(:,d+j),1,samplesize)'.*(zdiffall<=0),2)./expectDen);
    end
end
if (sum(isnan(scoreValall),'all')==0)
    checkError(6,1) = 'Score has NaN:';
    checkError(6,2) = 0;
else
    checkError(6,1) = 'Score has NaN:';
    checkError(6,2) = 1;
end
if (checkError(3,2)=='1' && checkError(1,2)=='0')
    checkError(7,1) = 'kbz or (2,2) is problematic';
    checkError(7,2) = 1;
end

scoreValall = (scoreValall);
end