function lambdaVec = lambdaD1EstHDM(parasLower,inputt,inputg,xall,zall,aall,deltaall,samplesize,goodI,d,bands)
estSample = length(inputg(:,1));
% read bandwidth and bandwidth adjustment
numeVec = zeros(samplesize,estSample,d,'single');
bandb = bands(1,1);
bandadjb = bands(1,3);
bandnh = bands(3,1);
bandadjnhn = bands(3,2);
bandadjnhd = bands(3,3);
% calculate all Z_i-estZ
aallmat = aall*ones(1,samplesize,'single');
zdiffall = zall*ones(1,samplesize,'single')-ones(samplesize,1,'single')*zall';% {i,j} = Z_i-Z_j
zdiffest = zall*ones(1,estSample,'single')-ones(samplesize,1,'single')*inputt';% {i,j} = Z_i-Z_j
kbz = kernel(zdiffest/bandb/bandadjb,1,'norm')/bandb/bandadjb;% {i,j} = Z_i-Z_j
% calculate \bb\trans\x
btransx = xall*[eye(d);parasLower]; % paras1 include upper d-by-d block.
% calculate all \bb\trans\x_i-\bb\trans\x_j (i row j column)
if (d>1)
    btxall = permute(repmat(inputg',[1 1 samplesize]),[3 2 1])-permute(repmat(btransx,[1 1 estSample]),[1 3 2]);% {i,j} = \bb\trans\x_j-\bb\trans\x_i
else
    btxall = repmat(btransx,[1 estSample])-repmat(inputg',[samplesize 1]);% {i,j} = \bb\trans\x_j-\bb\trans\x_i
end

% calculate nonparametric \lambda and \lambda\prime
if (d>1)
    kerD = 3;
else
    kerD = 1;
end
khbetax = kernel(btxall/bandnh/bandadjnhd,kerD,'norm')/(bandnh*bandadjnhd)^d;
denoVec = ((zdiffall<=0).*aallmat')*khbetax;
denoVec = denoVec(goodI,:);
% calculate kernel values of K\prime(.)
khprimevec = kernel2(btxall/bandnh/bandadjnhn,kerD,d,'quar')/(bandnh*bandadjnhn)^(d+1);
for t = 1:d
    numeVec(:,:,t) = ((zdiffall>=0).*aallmat)*squeeze(khprimevec(:,:,t));
end

lambdaVec = zeros(estSample,d,'single');
% khbetax = kernel(btxall/bandnh/bandadjnhd,kerD,'norm')/bandnh/bandadjnhd;
for t = 1:d
    lambdaVec(:,t) = sum(-kbz(goodI,:).*repmat(aall(goodI).*deltaall(goodI),1,estSample).*squeeze(khprimevec(goodI,:,t))./denoVec,1)'...
        +sum(kbz(goodI,:).*repmat(aall(goodI).*deltaall(goodI),1,estSample).*khbetax(goodI,:).*squeeze(numeVec(goodI,:,t))./(denoVec.^2),1)';
end

end